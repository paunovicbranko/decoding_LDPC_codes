The task of this project was to compare the performance of four different LDPC (Low-Density Parity-Check)
decoding algorithms: bit-flipping, Gallager-B, Gradient Descent Bit-flipping with Momentum, and without Momentum.

The first step involved generating codewords and passing them through a channel with errors. The errors were
implemented by inverting bits based on the random() function from the Python random library, with six different
error probabilities. This process was repeated 100,000 times.

We loaded the project's given parity-check matrix using the loadtxt() function from the numpy library.

• Bit-flipping:
This decoder was implemented by calculating all syndromes for each information bit and performing a majority
decision. Majority decision was implemented using the helper function decisionMaking(). If the result of the majority
decision was -1, the sign of the observed information bit was changed. This process was repeated for a predefined
maximum number of iterations or until the received sequence converged to the codeword. This algorithm is the simplest,
and therefore, it is expected to have the worst performance among all.

• Gradient Descent Bit-flipping with Momentum:
This algorithm represents an improved version of the bit-flipping algorithm. It is based on calculating the local
energy for each bit and using momentum to accelerate the optimization process.The threshold energy was chosen as the
minimum energy among all calculated energies for the bits. We compared the energy of each bit with the threshold energy,
and if it was less than or equal, and the probability of bit-flipping was less than 0.9, we inverted that bit.
When a bit is inverted, the ln variable, representing the momentum duration, is reset to zero. This process was repeated
for a predefined maximum number of iterations or until the value of all syndromes became 1. The values used as momentum
were obtained in the project settings.

• Gradient Descent Bit-flipping without Momentum:
This algorithm represents a simplified version of the previous algorithm. The procedure is almost the same, with the
difference that momentum is not used in the energy calculation. Therefore, variables L, ln, and rho were removed from
the code. Considering the removal of momentum, we expect this algorithm to have worse performance than the previous one.

• Gallager-B:
This algorithm represents a simple version of the message-passing algorithm, where nodes exchange one-bit messages.
Variable nodes perform a majority decision based on the received messages, except for the message received from the node
to which the message is being sent. Check nodes multiply the received messages, except for the message from the variable
node to which the message is being sent. The final message is obtained as the majority decision of each column of the
lambda_msgs matrix individually.
