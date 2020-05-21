# matrixprofile_js
An implementation of the UCR Matrix Profile and related work in JavaScript.

The Scalable Time series Anytime Matrix Profile (STAMP) algorithm:
----------------------------------------------
Takes as arguments two sequences (or the same sequence, twice, in the case of a self-join) and a subsequence length.
Returns the matrix profile and matrix profile index for the given sequence(s) and subsequence length.
The largest value(s) in the matrix profile corresponds to the most discordant pattern(s) (n discords).
The smallest values (in pairs or larger groups) in the matrix profile corresponds to the most similar patterns (n motifs).
There are other uses for the matrix profile as well (segmentation, measuring complexity, etc. etc. I really recommend following the work of the UCR team).

Mueen's Algorithm for Similarity Search (MASS):
----------------------------------------------
Takes as arguments a sequence and a query sequence.
Returns a distance profile which effectively represents a similarity search for the query over the given sequence.
Lower values will represent more similar patterns at that index

I'll also be adding the motif, discord, and segmentation methods from the paper shortly.

It is best to refer to the UCR Matrix Profile page for more details https://www.cs.ucr.edu/~eamonn/MatrixProfile.html

ComplexArray and all FFT JavaScript taken (literally copy and pasted) from https://github.com/dntj/jsfft/ which is under an MIT license.

I noted the following issues when testing and comparing with Numpy's FFT and introduced the following functions to correct them:
| Issue noticed | Function introduced |
|---------------|---------------------|
|Imaginary components of the FFT and invFFT appeared to have their sign inverted|     signCorrect     |
|output of the InvFFT appears not to be divided by the length of the resulting series as it should be|     invNormCorrect  |

The implementation of:
STAMP : Scalable Time series Anytime Matrix Profile
MASS : Mueen's algorithm for similarity search
are taken directly from the Matrix Profile I paper https://www.cs.ucr.edu/~eamonn/PID4481997_extend_Matrix%20Profile_I.pdf
All credit due to Eamonn Keogh, Abduallah Mueen, and their respective teams.

In validating and testing this implementation, I relied heavily on the Python matrixprofile library https://github.com/matrix-profile-foundation/matrixprofile as well as the dummy page available here https://ui.matrixprofile.org/ for reference and validation of results.
