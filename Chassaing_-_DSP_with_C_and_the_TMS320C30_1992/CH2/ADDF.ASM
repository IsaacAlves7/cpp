;ADDF.ASM - ADD 5 FLOATING-POINT VALUES
          .TITLE  "ADDF.ASM"     ;PROGRAM TITLE
          .GLOBAL BEGIN          ;REF/DEF SYMBOLS
VAL_ADDR  .WORD   VALUES         ;SPACES FOR VALUES      
VALUES    .FLOAT  0,1,2,3,4      ;5 VALUES TO BE ADDED
BEGIN     LDI     @VAL_ADDR,AR0  ;LOAD ADDR OF VALUES->AR0
          LDF     0.0,R0         ;INIT R0=0
          CALL    ADDF           ;CALL SUBROUTINE ADDF
WAIT      BR      WAIT           ;WAIT
;SUBROUTINE ADDF                 
ADDF      RPTS    4              ;NEXT INSTRUCTION 5 TIMES
          ADDF    *AR0++,R0      ;ACCUMULATE INTO R0
          RETS                   ;RETURN FROM SUBROUTINE
          .END                   ;END

