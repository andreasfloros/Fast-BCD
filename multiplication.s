              ;       BCDBIGMUL
              ;       Output: signed BCD multiplication memR2 = memR1 * memR0
              ;       R0 and R1 are assumed to be pointing to valid BCD numbers
              ;       of length R3 words (This subroutine works for R3 = 1, 2, 4, 8)

              ;       The subroutine works as follows:
              ;       We get the Most Significant Digits of memR0 and memR1 and check if they are negative.
              ;       We convert negative values to positive values via 10s complement inversion.
              ;       After this is done the problem is reduced to unsigned multiplication.
              ;       After unsigned multiplication is carried out we check if we need to negate/invert the result.
              ;       Inversion is required only when the two input numbers have different signs.
              ;       We implement this condition with an EOR of the MSDs.
              ;       If this condition is true we invert/negate otherwise we return the result as is.

              ;       This is the general approach that is taken.
              ;       In terms of optimizations, we focus on the unsigned multiplication algorithm.
              ;       This part of the code determines the overall efficiency.

              ;       Unsigned multiplication can be approached through various ways:
              ;       Classic long multiplication O(n^2)
              ;       Karatsuba multiplication ~O(n^1.58)
              ;       Toom Cook multiplication (Complexity depends on parameter choice)
              ;       FFT methods O(nlognloglogn)

              ;       Note that the above is slightly misleading.
              ;       64 digit multiplication (K = 8) isn't large enough for these algorithms to be practical.
              ;       For example, (in general) Karatsuba becomes practical past ~500 bit inputs. 64 digits are only ~216 bits long.
              ;       Furthermore, big O complexity ignores constant multiples and lower order terms (since it's asymptotic).

              ;       Nevertheless, variations of these algorithms may be applied to yield decent results.
              ;       We choose to implement Karatsuba multiplication as it is the simplest sub quadratic algorithm listed.
              ;       It is the simplest in the sense that it introduces the least amount of lower order terms (very few additions etc compared to the other algorithms listed).

              ;       We choose to implement this algorithm in stages avoiding recursion all together.
              ;       This is done because recursion as a concept is very inefficient when it comes to memory inputs and outputs.
              ;       If we were to implement it we would have to store and load hundereds of words every time.
              ;       This number would increase significantly for larger inputs. (Since the function calls would grow 3^log2(n) ~ n^1.58)

              ;       We split the unsigned multiplication in four layers:
              ;       The first, second and third layers are essentially indentical.
              ;       They are KaratsubaMUL, KaratsubaMUL1 and KaratsubaMUL2 in the code.
              ;       They simulate recursion.
              ;       The advantage of this approach is that unlike recursion, no memory fields are stored and retrieved between calls.
              ;       This significatnly boosts the efficiency of the algorithm. (If we wanted to extend the algorithm to larger numbers we could add more layers)
              ;       The fourth layer performs unsigned multiplication on single word memory inputs (8 digits).
              ;       This is usgnmul in the code.

              ;       Every time we call our unsigned multiplication algorithm we check if the input is one word long.
              ;       If it is we go to usgnmul, else we move through the first three layers until we reach a word multiplication.
              ;       We choose to do long multiplication once we reach 8 digits.
              ;       The reason for this choice is because if we chose to switch algorithms at say 16 digits we would get a worse performance. (Seen through experimentation)
              ;       On the other hand, we can see that 8 digits is the optimal length to do long multiplication by observing the following:
              ;       Without karatsuba (just long multiplication) the performance for K = 8 is ~600k cycles
              ;       Applying karatsuba once gives ~315k cycles
              ;       Applying karatsuba twice gives ~218k cycles
              ;       Applying karatsuba three times gives ~201k cycles
              ;       We see that while the performance is increasing, it is doing so at a decreasing rate.
              ;       This means that we are nearing a limit after which applying karatsuba will be inefficient.
              ;       Seeing as the performance gain from the last two stages is only 17k we expect minimal improvements from further stages.
              ;       If anything, they might prove to be less efficient given the extra additions, shifts etc
              ;       Having the base case be a single word multiplication also has another advantage, inputs can now be contained within a register.
              ;       This allows for some small optimizations such as removing loops that were previously necessary to get the next word of an input.

              ;       To further optimize our algorithm we unroll all loops and try to minimize branches.
              ;       This technique significantly improves the performance of the algorithm and it is primarily used in subroutines digitmul and addtor2.

              ;       The resulting algorithm is a hybrid long multiplication - karatsuba algorithm.

              ;       Note: While the above analysis is true for normal binary multiplication, BCD multiplication adds a considerable overhead, making a full karatsuba approach the optimal strategy.
              ;       Todo: Optimise the leaf operations with Karatsuba as well.

BCDBIGMUL     STMFD   SP!, {R0,R1,R3-R6,R10-R11,LR}
              MOV     R6, #0
              LSL     R3, R3, #2
              SUB     R3, R3, #4
              LDR     R4, [R0,R3]
              LDR     R5, [R1,R3]
              LSR     R4, R4, #28
              LSR     R5, R5, #28
              ADD     R3, R3, #4
              CMP     R4, #5
              ADC     R4, R6, #0
              BLHS    invertr0
              CMP     R5, #5
              ADC     R5, R6, #0
              BLHS    invertr1
              BL      karatsubaMUL
              EORS    R4, R4, R5
              BLNE    invertr2
              LDMFD   SP!, {R0,R1,R3-R6,R10-R11,PC}

karatsubaMUL  STMFD   SP!, {R0-R8,R10-R11,LR}
              BL      initr2
              CMP     R3, #4
              BHI     reduce
              BL      usgnmul
              LDMFD   SP!, {R0-R8,R10-R11,PC}

reduce        BL      clearfields0
              MOV     R8, R2
              LSR     R3, R3, #1
              BL      loadhalf0
              LDR     R0, =r0bottom0
              LDR     R1, =r1bottom0
              LDR     R2, =low0
              BL      karatsubaMUL1
              LDR     R0, =r0top0
              LDR     R1, =r1top0
              LDR     R2, =high0
              BL      karatsubaMUL1
              LDR     R0, =r1bottom0
              LDR     R2, =sum10
              BL      BCDADD
              LDR     R2, [R2,R3]
              ADD     R11, R6, R2
              LDR     R0, =r0bottom0
              LDR     R1, =r0top0
              LDR     R2, =sum00
              BL      BCDADD
              LDR     R2, [R2,R3]
              ADD     R10, R6, R2
              LDR     R0, =sum00
              LDR     R1, =sum10
              LDR     R2, =mid0
              BL      karatsubaMUL1
              ORRS    R7, R10, R11
              BLNE    carrymul
              LDR     R0, =low0
              LDR     R1, =high0
              LDR     R2, =lowandhigh0
              BL      BCDADD
              LSL     R3, R3, #1
              BL      invertr2
              LSR     R3, R3, #1
              LDR     R2, =mid0
              LDR     R0, =mid0
              LDR     R1, =lowandhigh0
              BL      BCDADD
              LDR     R11, =mid0
              LSL     R5, R3, #3
              LSL     R3, R3, #1
              BL      shift
              LSR     R3, R3, #1
              LDR     R1, =low0
              MOV     R2, R8
              BL      BCDADD
              MOV     R1, R2
              LDR     R0, =high0
              MOV     R11, R0
              LSL     R5, R3, #4
              LSL     R3, R3, #1
              BL      shift
              LSR     R3, R3, #1
              BL      BCDADD
              LDMFD   SP!, {R0-R8,R10-R11,PC}

karatsubaMUL1 STMFD   SP!, {R0-R8,R10-R11,LR}
              BL      initr2
              CMP     R3, #4
              BHI     reduce1
              BL      usgnmul
              LDMFD   SP!, {R0-R8,R10-R11,PC}

reduce1       BL      clearfields1
              MOV     R8, R2
              LSR     R3, R3, #1
              BL      loadhalf1
              LDR     R0, =r0bottom1
              LDR     R1, =r1bottom1
              LDR     R2, =low1
              BL      karatsubaMUL2
              LDR     R0, =r0top1
              LDR     R1, =r1top1
              LDR     R2, =high1
              BL      karatsubaMUL2
              LDR     R0, =r1bottom1
              LDR     R2, =sum11
              BL      BCDADD
              LDR     R2, [R2,R3]
              ADD     R11, R6, R2
              LDR     R0, =r0bottom1
              LDR     R1, =r0top1
              LDR     R2, =sum01
              BL      BCDADD
              LDR     R2, [R2,R3]
              ADD     R10, R6, R2
              LDR     R0, =sum01
              LDR     R1, =sum11
              LDR     R2, =mid1
              BL      karatsubaMUL2
              ORRS    R7, R10, R11
              BLNE    carrymul
              LDR     R0, =low1
              LDR     R1, =high1
              LDR     R2, =lowandhigh1
              BL      BCDADD
              LSL     R3, R3, #1
              BL      invertr2
              LSR     R3, R3, #1
              LDR     R2, =mid1
              LDR     R0, =mid1
              LDR     R1, =lowandhigh1
              BL      BCDADD
              LDR     R11, =mid1
              LSL     R5, R3, #3
              LSL     R3, R3, #1
              BL      shift
              LSR     R3, R3, #1
              LDR     R1, =low1
              MOV     R2, R8
              BL      BCDADD
              MOV     R1, R2
              LDR     R0, =high1
              MOV     R11, R0
              LSL     R5, R3, #4
              LSL     R3, R3, #1
              BL      shift
              LSR     R3, R3, #1
              BL      BCDADD
              LDMFD   SP!, {R0-R8,R10-R11,PC}

karatsubaMUL2 STMFD   SP!, {R0-R8,R10-R11,LR}
              BL      initr2
              CMP     R3, #4
              BHI     reduce2
              BL      usgnmul
              LDMFD   SP!, {R0-R8,R10-R11,PC}

reduce2       BL      clearfields2
              MOV     R8, R2
              LSR     R3, R3, #1
              BL      loadhalf2
              LDR     R0, =r0bottom2
              LDR     R1, =r1bottom2
              LDR     R2, =low2
              BL      karatsubaMUL2
              LDR     R0, =r0top2
              LDR     R1, =r1top2
              LDR     R2, =high2
              BL      karatsubaMUL2
              LDR     R0, =r1bottom2
              LDR     R2, =sum12
              BL      BCDADD
              LDR     R2, [R2,R3]
              ADD     R11, R6, R2
              LDR     R0, =r0bottom2
              LDR     R1, =r0top2
              LDR     R2, =sum02
              BL      BCDADD
              LDR     R2, [R2,R3]
              ADD     R10, R6, R2
              LDR     R0, =sum02
              LDR     R1, =sum12
              LDR     R2, =mid2
              BL      karatsubaMUL2
              ORRS    R7, R10, R11
              BLNE    carrymul
              LDR     R0, =low2
              LDR     R1, =high2
              LDR     R2, =lowandhigh2
              BL      BCDADD
              LSL     R3, R3, #1
              BL      invertr2
              LSR     R3, R3, #1
              LDR     R2, =mid2
              LDR     R0, =mid2
              LDR     R1, =lowandhigh2
              BL      BCDADD
              LDR     R11, =mid2
              LSL     R5, R3, #3
              LSL     R3, R3, #1
              BL      shift
              LSR     R3, R3, #1
              LDR     R1, =low2
              MOV     R2, R8
              BL      BCDADD
              MOV     R1, R2
              LDR     R0, =high2
              MOV     R11, R0
              LSL     R5, R3, #4
              LSL     R3, R3, #1
              BL      shift
              LSR     R3, R3, #1
              BL      BCDADD
              LDMFD   SP!, {R0-R8,R10-R11,PC}

carrymul      STMFD   SP!, {R0-R3,R6,LR}
              CMP     R10, #0
              STMFD   SP!, {R0,R2,R3}
              ADDNE   R2, R2, R3
              STR     R6, [R1,R3]
              SUBNE   R3, R3, R3, LSR #2
              MOVNE   R0, R2
              BLNE    BCDADD
              LDMFD   SP!, {R0,R2,R3}
              CMP     R11, #0
              STMFD   SP!, {R2,R3}
              ADDNE   R2, R2, R3
              STR     R6, [R0,R3]
              SUBNE   R3, R3, R3, LSR #2
              MOVNE   R1, R2
              BLNE    BCDADD
              LDMFD   SP!, {R2,R3}
              TST     R10, R11
              BNE     finaladd
              LDMFD   SP!, {R0-R3,R6,PC}

finaladd      SUBS    R10, R10, #0
              ADD     R2, R2, R3, LSL #1
              LDR     R6, [R2]
              BL      addc
              STR     R6, [R2]
              STRHS   R10, [R2,#4]
              LDMFD   SP!, {R0-R3,R6,PC}

clearfields0  STMFD   SP!, {R0-R8,R10,LR}
              MOV     R0, #0
              MOV     R2, #0
              MOV     R3, #0
              MOV     R4, #0
              MOV     R5, #0
              MOV     R7, #0
              MOV     R8, #0
              MOV     R10, #0
              LDR     R1, =r0top0
              B       clearfl

clearfields1  STMFD   SP!, {R0-R8,R10,LR}
              MOV     R0, #0
              MOV     R2, #0
              MOV     R3, #0
              MOV     R4, #0
              MOV     R5, #0
              MOV     R7, #0
              MOV     R8, #0
              MOV     R10, #0
              LDR     R1, =r0top1
              B       clearfl

clearfields2  STMFD   SP!, {R0-R8,R10,LR}
              MOV     R0, #0
              MOV     R2, #0
              MOV     R3, #0
              MOV     R4, #0
              MOV     R5, #0
              MOV     R7, #0
              MOV     R8, #0
              MOV     R10, #0
              LDR     R1, =r0top2
clearfl       STMIA   R1!, {R0,R2,R3,R4,R5,R7,R8,R10}
              ADD     R6, R6, #1
              CMP     R6, #20
              BNE     clearfl
              LDMFD   SP!, {R0-R8,R10,PC}

loadhalf0     STMFD   SP!, {R0,R1,R4,R5,LR}
              MOV     R4, #0
halflbottom0  LDR     R5, [R0,R4]
              STR     R5, [R4,#r0bottom0]
              LDR     R5, [R1,R4]
              STR     R5, [R4,#r1bottom0]
              ADD     R4, R4, #4
              CMP     R4, R3
              BNE     halflbottom0
              MOV     R4, #0
              ADD     R0, R0, R3
              ADD     R1, R1, R3
halfltop0     LDR     R5, [R0,R4]
              STR     R5, [R4,#r0top0]
              LDR     R5, [R1,R4]
              STR     R5, [R4,#r1top0]
              ADD     R4, R4, #4
              CMP     R4, R3
              BNE     halfltop0
              LDMFD   SP!, {R0,R1,R4,R5,PC}

loadhalf1     STMFD   SP!, {R0,R1,R4,R5,R6,R7,LR}
              MOV     R4, #0
              LDR     R6, =r0bottom1
              LDR     R7, =r1bottom1
halflbottom1  LDR     R5, [R0,R4]
              STR     R5, [R4,R6]
              LDR     R5, [R1,R4]
              STR     R5, [R4,R7]
              ADD     R4, R4, #4
              CMP     R4, R3
              BNE     halflbottom1
              MOV     R4, #0
              ADD     R0, R0, R3
              ADD     R1, R1, R3
              LDR     R6, =r0top1
              LDR     R7, =r1top1
halfltop1     LDR     R5, [R0,R4]
              STR     R5, [R4,R6]
              LDR     R5, [R1,R4]
              STR     R5, [R4,R7]
              ADD     R4, R4, #4
              CMP     R4, R3
              BNE     halfltop1
              LDMFD   SP!, {R0,R1,R4,R5,R6,R7,PC}

loadhalf2     STMFD   SP!, {R0,R1,R4,R5,R6,R7,LR}
              MOV     R4, #0
              LDR     R6, =r0bottom2
              LDR     R7, =r1bottom2
halflbottom2  LDR     R5, [R0,R4]
              STR     R5, [R4,R6]
              LDR     R5, [R1,R4]
              STR     R5, [R4,R7]
              ADD     R4, R4, #4
              CMP     R4, R3
              BNE     halflbottom2
              MOV     R4, #0
              LDR     R6, =r0top2
              LDR     R7, =r1top2
              ADD     R0, R0, R3
              ADD     R1, R1, R3
halfltop2     LDR     R5, [R0,R4]
              STR     R5, [R4,R6]
              LDR     R5, [R1,R4]
              STR     R5, [R4,R7]
              ADD     R4, R4, #4
              CMP     R4, R3
              BNE     halfltop2
              LDMFD   SP!, {R0,R1,R4,R5,R6,R7,PC}

initr2        STMFD   SP!, {R0,LR}
              MOV     R0, #0
init_loop     STR     R6, [R2,R0]
              ADD     R0, R0, #4
              CMP     R0, R3, LSL #1
              BNE     init_loop
              LDMFD   SP!, {R0,PC}

invertr0      STMFD   SP!, {R4-R7,LR}
              MOV     R4, #0
              LDR     R7, =inv_const
              LDR     R5, [R4,R7]
r0loop        LDR     R6, [R0,R4]
              SUB     R6, R5, R6
              BL      addc
              STR     R6, [R4,#r0field]
              ADD     R4, R4, #4
              TEQ     R4, R3
              BNE     r0loop
              LDR     R0, =r0field
              LDMFD   SP!, {R4-R7,PC}

invertr1      STMFD   SP!, {R4-R7,LR}
              MOV     R4, #0
              LDR     R7, =inv_const
              LDR     R5, [R4,R7]
r1loop        LDR     R6, [R1,R4]
              SUB     R6, R5, R6
              BL      addc
              STR     R6, [R4,#r1field]
              ADD     R4, R4, #4
              TEQ     R4, R3
              BNE     r1loop
              LDR     R1, =r1field
              LDMFD   SP!, {R4-R7,PC}

invertr2      STMFD   SP!, {R3-R7,LR}
              SUBS    R4, R6, #0
              LSL     R3, R3, #1
              LDR     R7, =inv_const
              LDR     R5, [R4,R7]
r2loop        LDR     R6, [R2,R4]
              SUB     R6, R5, R6
              BL      addc
              STR     R6, [R2,R4]
              ADD     R4, R4, #4
              TEQ     R4, R3
              BNE     r2loop
              LDMFD   SP!, {R3-R7,PC}

addc          STMFD   SP!, {R2,R5,R10,LR}
              MOV     R10, #0xF
              AND     R5, R10, R6
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              MOV     R2, R5
              AND     R5, R10, R6, LSR #4
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R2, R2, R5, LSL #4
              AND     R5, R10, R6, LSR #8
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R2, R2, R5, LSL #8
              AND     R5, R10, R6, LSR #12
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R2, R2, R5, LSL #12
              AND     R5, R10, R6, LSR #16
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R2, R2, R5, LSL #16
              AND     R5, R10, R6, LSR #20
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R2, R2, R5, LSL #20
              AND     R5, R10, R6, LSR #24
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R2, R2, R5, LSL #24
              AND     R5, R10, R6, LSR #28
              ADC     R5, R5, #0
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R2, R2, R5, LSL #28
              MOV     R6, R2
              LDMFD   SP!, {R2,R5,R10,PC}

usgnmul       STMFD   SP!, {R4-R8,R11,LR}
              BL      clearsum
              LDR     R11, =sum
              MOV     R5, #0
              MOV     R6, #0xF
r0words       LDR     R7, [R0]
              AND     R8, R6, R7
              BL      digitmul
              BL      shift1
              BL      addtor2
              ADD     R5, R5, #4
              AND     R8, R6, R7, LSR #4
              BL      digitmul
              BL      shift1
              BL      addtor2
              ADD     R5, R5, #4
              AND     R8, R6, R7, LSR #8
              BL      digitmul
              BL      shift1
              BL      addtor2
              ADD     R5, R5, #4
              AND     R8, R6, R7, LSR #12
              BL      digitmul
              BL      shift1
              BL      addtor2
              ADD     R5, R5, #4
              AND     R8, R6, R7, LSR #16
              BL      digitmul
              BL      shift1
              BL      addtor2
              ADD     R5, R5, #4
              AND     R8, R6, R7, LSR #20
              BL      digitmul
              BL      shift1
              BL      addtor2
              ADD     R5, R5, #4
              AND     R8, R6, R7, LSR #24
              BL      digitmul
              BL      shift1
              BL      addtor2
              ADD     R5, R5, #4
              AND     R8, R6, R7, LSR #28
              BL      digitmul
              BL      shift1
              BL      addtor2
              LDMFD   SP!, {R4-R8,R11,PC}

clearsum      STMFD   SP!, {R0-R2,LR}
              LDR     R0, =sum
              MOV     R1, #0
              MOV     R2, #0
clearl        STR     R2, [R0,R1]
              ADD     R1, R1, #4
              CMP     R1, #64
              BNE     clearl
              LDMFD   SP!, {R0-R2,PC}

digitmul      STMFD   SP!, {R2,R4-R7,R10,LR}
              LDR     R2, =mul0
              MOV     R4, #0
r1words       LDR     R5, [R1]
              AND     R6, R5, #0xF
              LSL     R7, R6, #3
              ADD     R7, R7, R6, LSL #2
              ADD     R7, R2, R7
              LDRB    R6, [R7,R8]
              AND     R7, R6, #0xF
              ADD     R7, R7, R4
              CMP     R7, #10
              SUBHS   R7, R7, #10
              ADDHS   R6, R6, #0x10
              LSR     R4, R6, #4
              MOV     R6, R7
              AND     R7, R5, #0xF0
              LSR     R7, R7, #4
              LSL     R10, R7, #3
              ADD     R10, R10, R7, LSL #2
              ADD     R10, R2, R10
              LDRB    R7, [R10,R8]
              AND     R10, R7, #0xF
              ADD     R10, R10, R4
              CMP     R10, #10
              SUBHS   R10, R10, #10
              ADDHS   R7, R7, #0x10
              LSR     R4, R7, #4
              ADD     R6, R6, R10, LSL #4
              AND     R7, R5, #0xF00
              LSR     R7, R7, #8
              LSL     R10, R7, #3
              ADD     R10, R10, R7, LSL #2
              ADD     R10, R2, R10
              LDRB    R7, [R10,R8]
              AND     R10, R7, #0xF
              ADD     R10, R10, R4
              CMP     R10, #10
              SUBHS   R10, R10, #10
              ADDHS   R7, R7, #0x10
              LSR     R4, R7, #4
              ADD     R6, R6, R10, LSL #8
              AND     R7, R5, #0xF000
              LSR     R7, R7, #12
              LSL     R10, R7, #3
              ADD     R10, R10, R7, LSL #2
              ADD     R10, R2, R10
              LDRB    R7, [R10,R8]
              AND     R10, R7, #0xF
              ADD     R10, R10, R4
              CMP     R10, #10
              SUBHS   R10, R10, #10
              ADDHS   R7, R7, #0x10
              LSR     R4, R7, #4
              ADD     R6, R6, R10, LSL #12
              AND     R7, R5, #0xF0000
              LSR     R7, R7, #16
              LSL     R10, R7, #3
              ADD     R10, R10, R7, LSL #2
              ADD     R10, R2, R10
              LDRB    R7, [R10,R8]
              AND     R10, R7, #0xF
              ADD     R10, R10, R4
              CMP     R10, #10
              SUBHS   R10, R10, #10
              ADDHS   R7, R7, #0x10
              LSR     R4, R7, #4
              ADD     R6, R6, R10, LSL #16
              AND     R7, R5, #0xF00000
              LSR     R7, R7, #20
              LSL     R10, R7, #3
              ADD     R10, R10, R7, LSL #2
              ADD     R10, R2, R10
              LDRB    R7, [R10,R8]
              AND     R10, R7, #0xF
              ADD     R10, R10, R4
              CMP     R10, #10
              SUBHS   R10, R10, #10
              ADDHS   R7, R7, #0x10
              LSR     R4, R7, #4
              ADD     R6, R6, R10, LSL #20
              AND     R7, R5, #0xF000000
              LSR     R7, R7, #24
              LSL     R10, R7, #3
              ADD     R10, R10, R7, LSL #2
              ADD     R10, R2, R10
              LDRB    R7, [R10,R8]
              AND     R10, R7, #0xF
              ADD     R10, R10, R4
              CMP     R10, #10
              SUBHS   R10, R10, #10
              ADDHS   R7, R7, #0x10
              LSR     R4, R7, #4
              ADD     R6, R6, R10, LSL #24
              AND     R7, R5, #0xF0000000
              LSR     R7, R5, #28
              LSL     R10, R7, #3
              ADD     R10, R10, R7, LSL #2
              ADD     R10, R2, R10
              LDRB    R7, [R10,R8]
              AND     R10, R7, #0xF
              ADD     R10, R10, R4
              CMP     R10, #10
              SUBHS   R10, R10, #10
              ADDHS   R7, R7, #0x10
              LSR     R4, R7, #4
              ADD     R6, R6, R10, LSL #28
              STR     R6, [R11]
              STR     R4, [R11,#4]
              LDMFD   SP!, {R2,R4-R7,R10,PC}

shift1        STMFD   SP!, {R0,R6,R3}
              ADD     R0, R3, #4
              MOV     R6, #0
              ADD     R3, R3, #2
init_loop1    STR     R6, [R11,R0]
              ADD     R0, R0, #4
              CMP     R0, R3, LSL #1
              BNE     init_loop1
              LDMFD   SP!, {R0,R6,R3}
shift         STMFD   SP!, {R0,R2,R3,R4,R5,R6,R7,LR}
              LSL     R3, R3, #1
nextshift     TEQ     R5, #0
              BEQ     escape
              MOV     R6, #0
              MOV     R0, #0
              MOV     R2, #0xF
sumwords      LDR     R4, [R0,R11]
              AND     R7, R2, R4, LSR #28
              LSL     R4, R4, #4
              ADD     R4, R4, R6
              STR     R4, [R0,R11]
              MOV     R6, R7
              ADD     R0 ,R0, #4
              CMP     R0, R3
              BNE     sumwords
              SUB     R5, R5, #4
              B       nextshift
escape        LDMFD   SP!, {R0,R2,R3,R4,R5,R6,R7,PC}

addtor2       STMFD   SP!, {R0,R2,R4,R5,R6,R7,R10,R11,LR}
              MOV     R10, #0xF
              ADDS    R10, R10, #0
              LDR     R4, [R11], #4
              LDR     R0, [R2], #4
              AND     R6, R10, R4
              AND     R5, R10, R0
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              MOV     R7, R5
              AND     R6, R10, R4, LSR #4
              AND     R5, R10, R0, LSR #4
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #4
              AND     R6, R10, R4, LSR #8
              AND     R5, R10, R0, LSR #8
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #8
              AND     R6, R10, R4, LSR #12
              AND     R5, R10, R0, LSR #12
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #12
              AND     R6, R10, R4, LSR #16
              AND     R5, R10, R0, LSR #16
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #16
              AND     R6, R10, R4, LSR #20
              AND     R5, R10, R0, LSR #20
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #20
              AND     R6, R10, R4, LSR #24
              AND     R5, R10, R0, LSR #24
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #24
              AND     R6, R10, R4, LSR #28
              AND     R5, R10, R0, LSR #28
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #28
              STR     R7, [R2,#-4]
              LDR     R4, [R11], #4
              LDR     R0, [R2], #4
              AND     R6, R10, R4
              AND     R5, R10, R0
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              MOV     R7, R5
              AND     R6, R10, R4, LSR #4
              AND     R5, R10, R0, LSR #4
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #4
              AND     R6, R10, R4, LSR #8
              AND     R5, R10, R0, LSR #8
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #8
              AND     R6, R10, R4, LSR #12
              AND     R5, R10, R0, LSR #12
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #12
              AND     R6, R10, R4, LSR #16
              AND     R5, R10, R0, LSR #16
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #16
              AND     R6, R10, R4, LSR #20
              AND     R5, R10, R0, LSR #20
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #20
              AND     R6, R10, R4, LSR #24
              AND     R5, R10, R0, LSR #24
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #24
              AND     R6, R10, R4, LSR #28
              AND     R5, R10, R0, LSR #28
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #28
              STR     R7, [R2,#-4]
              LDMFD   SP!, {R0,R2,R4,R5,R6,R7,R10,R11,PC}

BCDADD        STMFD   SP!, {R0,R1,R3,R4,R5,R6,R7,R8,R10,R11,LR}
              LSL     R3, R3, #2
              MOV     R10, #0xF
              LSRS    R8, R10, #5
next          LDR     R4, [R0], #4
              LDR     R11, [R1], #4
              AND     R6, R10, R4
              AND     R5, R10, R11
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              MOV     R7, R5
              AND     R6, R10, R4, LSR #4
              AND     R5, R10, R11, LSR #4
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #4
              AND     R6, R10, R4, LSR #8
              AND     R5, R10, R11, LSR #8
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #8
              AND     R6, R10, R4, LSR #12
              AND     R5, R10, R11, LSR #12
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #12
              AND     R6, R10, R4, LSR #16
              AND     R5, R10, R11, LSR #16
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #16
              AND     R6, R10, R4, LSR #20
              AND     R5, R10, R11, LSR #20
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #20
              AND     R6, R10, R4, LSR #24
              AND     R5, R10, R11, LSR #24
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #24
              AND     R6, R10, R4, LSR #28
              AND     R5, R10, R11, LSR #28
              ADC     R5, R6, R5
              CMP     R5, #10
              SUBHS   R5, R5, #10
              ADD     R7, R7, R5, LSL #28
              STR     R7, [R2,R8]
              ADD     R8, R8, #4
              TEQ     R8, R3
              BNE     next
              LDMFD   SP!, {R0,R1,R3,R4,R5,R6,R7,R8,R10,R11,PC}

r0field       FILL    64
r1field       FILL    64
r0top0        FILL    64
r0bottom0     FILL    64
r1top0        FILL    64
r1bottom0     FILL    64
lowandhigh0   FILL    64
low0          FILL    64
high0         FILL    64
mid0          FILL    64
sum00         FILL    64
sum10         FILL    64
r0top1        FILL    64
r0bottom1     FILL    64
r1top1        FILL    64
r1bottom1     FILL    64
lowandhigh1   FILL    64
low1          FILL    64
high1         FILL    64
mid1          FILL    64
sum01         FILL    64
sum11         FILL    64
r0top2        FILL    64
r0bottom2     FILL    64
r1top2        FILL    64
r1bottom2     FILL    64
lowandhigh2   FILL    64
low2          FILL    64
high2         FILL    64
mid2          FILL    64
sum02         FILL    64
sum12         FILL    64
sum           FILL    64
inv_const     DCD     0x99999999
mul0          DCB     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
mul1          DCB     0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0
mul2          DCB     0, 2, 4, 6, 8, 0x10, 0x12, 0x14, 0x16, 0x18, 0, 0
mul3          DCB     0, 3, 6, 9, 0x12, 0x15, 0x18, 0x21, 0x24, 0x27, 0, 0
mul4          DCB     0, 4, 8, 0x12, 0x16, 0x20, 0x24, 0x28, 0x32, 0x36, 0, 0
mul5          DCB     0, 5, 0x10, 0x15, 0x20, 0x25, 0x30, 0x35, 0x40, 0x45, 0, 0
mul6          DCB     0, 6, 0x12, 0x18, 0x24, 0x30, 0x36, 0x42, 0x48, 0x54, 0, 0
mul7          DCB     0, 7, 0x14, 0x21, 0x28, 0x35, 0x42, 0x49, 0x56, 0x63, 0, 0
mul8          DCB     0, 8, 0x16, 0x24, 0x32, 0x40, 0x48, 0x56, 0x64, 0x72, 0, 0
mul9          DCB     0, 9, 0x18, 0x27, 0x36, 0x45, 0x54, 0x63, 0x72, 0x81, 0, 0