# Fast-BCD

This is an optimised library for BCD operations implemented in ARM Assembly.

## Optimised addition

BCD addition is implemented in the ```addition.s``` file.

This implementation focuses on optimising the space complexity.

__BCD8ADD__: Assume 8 digit BCD numbers x and y in R0 and R1 respectively.

Then at the end of the subroutine:

R2 will contain z = x + y (BCD addition).

R3 will contain the signed overflow (0 or 1).

R4 will contain the unsigned overflow (0 or 1).

__BCD8ADDMEM__: Same operation as above except R0, R1 and R3 are pointers to the values x, y and z respectively.

## Optimised multiplication

BCD multiplication is implemented in the ```multiplication.s``` file.

This implementation focuses on optimising the time complexity.

__BCDBIGMUL__: Define R3 = K = 1, 2 ,4 ,8 and let R0 and R1 point to BCD numbers x and y respectively that are 8K digits long. We use R2 as a pointer to the output of the multiplication.

Then at the end of the subroutine:

R2 will be pointing to the value z = x * y (BCD multiplication).

Note that this will generate a 16K digit number.

## Tests

Tests can be found in the ```tests``` folder for both ```addition.s``` and ```multiplication.s```.
