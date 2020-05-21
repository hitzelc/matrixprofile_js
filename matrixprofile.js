/*
ComplexArray and all FFT JavaScript taken from https://github.com/dntj/jsfft/
with minor modifications.
I noted the following issues when testing and comparing with Numpy's FFT:
--FUNCTION INTRODUCED--:--ISSUE CORRECTED BY FUNCTION---------------------------------
* signCorrect  : Imaginary components of the FFT and invFFT appeared to have their sign inverted (hence signCorrect)
* invNormCorrect : output of the InvFFT appears not to be divided by the length of the resulting series as it should be

The implementation of:
STAMP : Scalable Time series Anytime Matrix Profile
MASS : Mueen's algorithm for similarity search
are taken from the Matrix Profile I paper https://www.cs.ucr.edu/~eamonn/PID4481997_extend_Matrix%20Profile_I.pdf
with heavy reliance on the fine Python matrixprofile library here https://github.com/matrix-profile-foundation/matrixprofile
for reference and validation of results (the moving mean and standard deviation is a near-direct translation)
*/
class baseComplexArray {
    constructor(other, arrayType = Float32Array) {
    if (other instanceof baseComplexArray) {
        // Copy constuctor.
        this.ArrayType = other.ArrayType;
        this.real = new this.ArrayType(other.real);
        this.imag = new this.ArrayType(other.imag);
    } else {
        this.ArrayType = arrayType;
        // other can be either an array or a number.
        this.real = new this.ArrayType(other);
        this.imag = new this.ArrayType(this.real.length);
    }

    this.length = this.real.length;
    }

    toString() {
    const components = [];

    this.forEach((value, i) => {
        components.push(
        `(${value.real.toFixed(2)}, ${value.imag.toFixed(2)})`
        );
    });

    return `[${components.join(', ')}]`;
    }

    forEach(iterator) {
    const n = this.length;
    // For gc efficiency, re-use a single object in the iterator.
    const value = Object.seal(Object.defineProperties({}, {
        real: {writable: true}, imag: {writable: true},
    }));

    for (let i = 0; i < n; i++) {
        value.real = this.real[i];
        value.imag = this.imag[i];
        iterator(value, i, n);
    }
    }

    // In-place mapper.
    map(mapper) {
    this.forEach((value, i, n) => {
        mapper(value, i, n);
        this.real[i] = value.real;
        this.imag[i] = value.imag;
    });

    return this;
    }

    conjugate() {
    return new ComplexArray(this).map((value) => {
        value.imag *= -1;
    });
    }

    magnitude() {
    const mags = new this.ArrayType(this.length);

    this.forEach((value, i) => {
        mags[i] = Math.sqrt(value.real*value.real + value.imag*value.imag);
    });

    return mags;
    }
}

// Math constants and functions we need.
const PI = Math.PI;
const SQRT1_2 = Math.SQRT1_2;

function FFT(input) {
    return ensureComplexArray(input).FFT();
}

function InvFFT(input) {
    return ensureComplexArray(input).InvFFT();
}

function frequencyMap(input, filterer) {
    return ensureComplexArray(input).frequencyMap(filterer);
}

class ComplexArray extends baseComplexArray {
    FFT() {
    return fft(this, false);
    }

    InvFFT() {
    return fft(this, true);
    }

    // Applies a frequency-space filter to input, and returns the real-space
    // filtered input.
    // filterer accepts freq, i, n and modifies freq.real and freq.imag.
    frequencyMap(filterer) {
    return this.FFT().map(filterer).InvFFT();
    }
}

function ensureComplexArray(input) {
    return input instanceof ComplexArray && input || new ComplexArray(input);
}

function fft(input, inverse) {
    const n = input.length;

    if (n & (n - 1)) {
    return FFT_Recursive(input, inverse);
    } else {
    return FFT_2_Iterative(input, inverse);
    }
}

function FFT_Recursive(input, inverse) {
    const n = input.length;

    if (n === 1) {
    return input;
    }

    const output = new ComplexArray(n, input.ArrayType);

    // Use the lowest odd factor, so we are able to use FFT_2_Iterative in the
    // recursive transforms optimally.
    const p = LowestOddFactor(n);
    const m = n / p;
    const normalisation = 1 / Math.sqrt(p);
    let recursive_result = new ComplexArray(m, input.ArrayType);

    // Loops go like O(n S p_i), where p_i are the prime factors of n.
    // for a power of a prime, p, this reduces to O(n p log_p n)
    for(let j = 0; j < p; j++) {
    for(let i = 0; i < m; i++) {
        recursive_result.real[i] = input.real[i * p + j];
        recursive_result.imag[i] = input.imag[i * p + j];
    }
    // Don't go deeper unless necessary to save allocs.
    if (m > 1) {
        recursive_result = fft(recursive_result, inverse);
    }

    const del_f_r = Math.cos(2*PI*j/n);
    const del_f_i = (inverse ? -1 : 1) * Math.sin(2*PI*j/n);
    let f_r = 1;
    let f_i = 0;

    for(let i = 0; i < n; i++) {
        const _real = recursive_result.real[i % m];
        const _imag = recursive_result.imag[i % m];

        output.real[i] += f_r * _real - f_i * _imag;
        output.imag[i] += f_r * _imag + f_i * _real;

        [f_r, f_i] = [
        f_r * del_f_r - f_i * del_f_i,
        f_i = f_r * del_f_i + f_i * del_f_r,
        ];
    }
    }
    /*
    for(let i=0;i<output.imag.length;i++){
    output.imag[i] *= -1;
    }*/
    
    // Copy back to input to match FFT_2_Iterative in-placeness
    /*
    if(inverse){
    for(let i = 0; i < n; i++) {
        output.real[i] = normalisation * output.real[i];
        output.imag[i] = normalisation * output.imag[i];
    }
    }*/
    /*
    if(inverse){
    for(let i = 0; i < output.length; i++) {
        output.real[i] = output.real[i]/output.length;
        output.imag[i] = output.imag[i]/output.length;
    }
    }*/

    return output;
}

function FFT_2_Iterative(input, inverse) {
    const n = input.length;

    const output = BitReverseComplexArray(input);
    const output_r = output.real;
    const output_i = output.imag;
    // Loops go like O(n log n):
    //   width ~ log n; i,j ~ n
    let width = 1;
    while (width < n) {
    const del_f_r = Math.cos(PI/width);
    const del_f_i = (inverse ? -1 : 1) * Math.sin(PI/width);
    for (let i = 0; i < n/(2*width); i++) {
        let f_r = 1;
        let f_i = 0;
        for (let j = 0; j < width; j++) {
        const l_index = 2*i*width + j;
        const r_index = l_index + width;

        const left_r = output_r[l_index];
        const left_i = output_i[l_index];
        const right_r = f_r * output_r[r_index] - f_i * output_i[r_index];
        const right_i = f_i * output_r[r_index] + f_r * output_i[r_index];

        output_r[l_index] = SQRT1_2 * (left_r + right_r);
        output_i[l_index] = SQRT1_2 * (left_i + right_i);
        output_r[r_index] = SQRT1_2 * (left_r - right_r);
        output_i[r_index] = SQRT1_2 * (left_i - right_i);

        [f_r, f_i] = [
            f_r * del_f_r - f_i * del_f_i,
            f_r * del_f_i + f_i * del_f_r,
        ];
        }
    }
    width <<= 1;
    }

    return output;
}

function BitReverseIndex(index, n) {
    let bitreversed_index = 0;

    while (n > 1) {
    bitreversed_index <<= 1;
    bitreversed_index += index & 1;
    index >>= 1;
    n >>= 1;
    }
    return bitreversed_index;
}

function BitReverseComplexArray(array) {
    const n = array.length;
    const flips = new Set();

    for(let i = 0; i < n; i++) {
    const r_i = BitReverseIndex(i, n);

    if(flips.has(i)){continue;}

    [array.real[i], array.real[r_i]] = [array.real[r_i], array.real[i]];
    [array.imag[i], array.imag[r_i]] = [array.imag[r_i], array.imag[i]];

    flips.add(r_i);
    }

    return array;
}

function LowestOddFactor(n) {
    const sqrt_n = Math.sqrt(n);
    let factor = 3;

    while(factor <= sqrt_n) {
    if (n % factor === 0) return factor;
    factor += 2;
    }
    return n;
}

function signCorrect(fftout){
    for(let i=0;i<fftout.imag.length;i++){
        fftout.imag[i] *= -1;
    }
    return fftout;
}

function invNormCorrect(ifftout){
    //console.log(ifftout.real.length);
    for(let i=0;i<ifftout.real.length;i++){
        ifftout.real[i] /= ifftout.real.length;
        ifftout.imag[i] /= ifftout.imag.length;
    }
    return ifftout;
}

function slidingDotProduct(q,t,runnum){
    const query = [...q];
    const timeseries = [...t];
    const m = query.length;
    const n = timeseries.length;
    query.reverse();
    //for(let i=0;i<n;i++){timeseries.push(0);}
    //for(let i=0;i<(2*n - m);i++){reversequery.push(0);}
    for(let i=0;i<(n - m);i++){query.push(0);}
    const reversequery_fourier = signCorrect(FFT(query));
    const timeseries_fourier = signCorrect(FFT(timeseries));
    const mult_query_timeseries = new ComplexArray(timeseries_fourier.length,timeseries_fourier.ArrayType);
    for(let i=0;i<timeseries_fourier.length;i++){
        mult_query_timeseries.real[i] = reversequery_fourier.real[i]*timeseries_fourier.real[i] - reversequery_fourier.imag[i]*timeseries_fourier.imag[i];
        mult_query_timeseries.imag[i] = reversequery_fourier.real[i]*timeseries_fourier.imag[i] + reversequery_fourier.imag[i]*timeseries_fourier.real[i];
    }
    //(reversequery_fourier.real[i] + reversequery_fourier.imag[i]) * (timeseries_fourier.real[i] + timeseries_fourier.imag[i])
    //= reversequery_fourier.real[i]*timeseries_fourier.real[i] + reversequery_fourier.real[i]*timeseries_fourier.imag[i] + reversequery_fourier.imag[i]*timeseries_fourier.real[i] + reversequery_fourier.imag[i]*timeseries_fourier.imag[i];
    //Re = rev_r*tim_r + rev_i*tim_i
    //Im = rev_i*tim_r + rev_r*tim_i
    const multinv = signCorrect(invNormCorrect(InvFFT(mult_query_timeseries)));
    const multinv_real = Array.from(multinv.real);
    //console.log(multinv_real);
    multinv_real.reverse();
    let firstelem = multinv_real.pop();
    multinv_real.unshift(firstelem);
    //if(runnum === 0){console.log(multinv_real);}
    return multinv_real.slice(m-1,n);
}

function apply_exclusion_zone(exclusion_zone, is_join, window_size, data_length, index, distance_profile){
    if(exclusion_zone > 0 && !is_join){
        ez_start = Math.max(0, index - exclusion_zone);
        ez_end = Math.min(data_length - window_size + 1, index + exclusion_zone);
        for(let i=0;i<distance_profile.length;i++){
            if(i >= ez_start && i < ez_end){
                distance_profile[i] = 1000000000000;
            }
        }
    }
    return distance_profile;
}

function distanceProfile(q,t,qt,Mu_t,Sd_t,mu_q,sd_q){
    const m = q.length;
    const n = t.length;
    const dist = [];
    for(let i=0;i<qt.length;i++){
        dist.push(Math.sqrt(2 * m * (1 - ((qt[i] - (m * Mu_t[i] * mu_q)) / (m* Sd_t[i]*sd_q)))));
    }
    return dist;
}

function add(a,b){return a+b;}

function mean(arr){
    const out = arr.reduce(add,0)/arr.length;
    return out;
}

function sdev(x){
    const avgx = mean(x);
    let xx = [];
    for(let i=0;i<x.length;i++){
        xx[i] = (x[i] - avgx)*(x[i] - avgx);
    }
    const sumxx = xx.reduce(add,0);
    const stdev = Math.sqrt( (sumxx / x.length) );
    return stdev;
}

function moving_mu_std(arr,window_size){
    //largely just a translation of https://github.com/matrix-profile-foundation/matrixprofile/blob/fd930ba54835bab07065449b3fdde2793ba395d4/matrixprofile/cycore.pyx
    const n = arr.length;
    const w = window_size;
    const profile_len = n - w + 1;
    const cumsum = new Array(n);
    const sq_cumsum = new Array(n);
    const sums = new Array(profile_len);
    const sq_sums = new Array(profile_len);
    const mu = new Array(profile_len);
    const sig_sq = new Array(profile_len);
    const sig = new Array(profile_len);
    cumsum[0] = arr[0];
    for(let i=1;i<n;i++){
        cumsum[i] = arr[i] + cumsum[i - 1];
    }
    sq_cumsum[0] = arr[0] * arr[0];
    for(let i=1;i<n;i++){
        sq_cumsum[i] = arr[i] * arr[i] + sq_cumsum[i - 1];
    }
    sums[0] = cumsum[w - 1];
    for(let i=0;i<n-w;i++){
        sums[i + 1] = cumsum[w + i] - cumsum[i];
    }
    sq_sums[0] = sq_cumsum[w - 1];
    for(let i=0;i<n-w;i++){
        sq_sums[i + 1] = sq_cumsum[w + i] - sq_cumsum[i];
    }
    for(let i=0;i<profile_len;i++){
        mu[i] = sums[i] / w;
    }
    for(let i=0;i<profile_len;i++){
        sig_sq[i] = (sq_sums[i] / w) - (mu[i] * mu[i]);
    }
    for(let i=0;i<profile_len;i++){
        sig[i] = Math.sqrt(sig_sq[i]);
    }
    return [mu, sig];
}

function MASS(query,timeseries,runnum){
    const qt = slidingDotProduct(query,timeseries,runnum);
    const mean_query = mean(query);
    const stdev_query = sdev(query);
    const movmustd_ts = moving_mu_std(timeseries,query.length);
    const mean_timeseries = movmustd_ts[0]; //should be length 756
    const stdev_timeseries = movmustd_ts[1];
    return distanceProfile(query,timeseries,qt,mean_timeseries,stdev_timeseries,mean_query,stdev_query);
}

function STAMP(t_a,t_b,subseq_length){
    const timeseries_a = Float32Array.from(t_a);
    const timeseries_b = Float32Array.from(t_b);
    const n_b = timeseries_b.length;
    const exclusion_zone = parseInt(Math.ceil(subseq_length / 2.0));
    const indexes = [];
    const matrix_profile = []; //inf all the way
    const matrix_profile_indexes = []; //0 all the way
    for(let i=0;i<(n_b - subseq_length);i++){
        matrix_profile.push(1000000000000);
        matrix_profile_indexes.push(0);
        indexes.push(i);
    }
    for(let i of indexes){
        const distprof = MASS(timeseries_b.slice(i,i+subseq_length),timeseries_a,i);
        const D = apply_exclusion_zone(exclusion_zone, false, subseq_length, timeseries_a.length, i, distprof);
        for(let d=0;d<D.length-1;d++){
            if(D[d] <= matrix_profile[d]){
                matrix_profile[d] = D[d];
                matrix_profile_indexes[d] = i;
            }
            else{continue;}
        }
    }
    return {'matrix_profile':matrix_profile,'matrix_profile_index':matrix_profile_indexes};
}
