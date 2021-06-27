           ; Perform BCD addition
           ; This implementation is optimised with respect to space and makes use of unconventional looping to achieve the best results
           
BCD8ADD    STMFD   SP!, {R5-R8,R10,LR}
           MOV     R2, #0
           ADDS    R8, R2, #0
           MOV     R10, #0xF
dg         AND     R6, R10, R0, LSR R8
           AND     R7, R10, R1, LSR R8
           ADC     R5, R6, R7
           CMP     R5, #10
           SUBHS   R5, R5, #10
           ADD     R2, R2, R5, LSL R8
           ADD     R8, R8, #4
           EORS    R3, R8, #32
           BNE     dg
           ADC     R4, R3, #0
           CMP     R6, #5
           ADC     R6, R3, #0
           CMP     R7, #5
           ADC     R7, R3, #0
           CMP     R5, #5
           ADC     R5, R3, #0
           CMP     R6, R7
           EOREQ   R3, R5, R6
           LDMFD   SP!, {R5-R8,R10,PC}
BCD8ADDMEM STMFD   SP!, {R0-R2,R5,LR}
           MOV     R5, R2
           LDR     R0, [R0]
           LDR     R1, [R1]
           BL      BCD8ADD
           STR     R2, [R5]
           LDMFD   SP!, {R0-R2,R5,PC}