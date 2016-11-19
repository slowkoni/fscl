%{

#include "msparser.h"
#include "ms-parser.h"
#ifdef __cplusplus
extern "C" {
#endif
%}
%option reentrant bison-bridge yylineno noyywrap
%x MS_BLOCK
%x EAT
%%
<EAT>{
^\/\/.*\n { BEGIN(MS_BLOCK); return MS_START; }
   ^.*\n { }
}
<MS_BLOCK>{
	\[		 { return OPEN_BRACKET;		}
	\]		 { return CLOSE_BRACKET;	}
	\(		 { return OPEN_PAREN;		}
	:		 { return COLON;		}
	\)		 { return CLOSE_PAREN;		}
	,		 { return COMMA;		}
	;		 { return SEMICOLON;		}
	^[01]+$		 { return HAPLOTYPE;		}
	[0-9]+		 { return INTEGER;		}
	[0-9]*(\.[0-9]+){0,1}(e[\+\-][0-9]+){0,1}	 { return DECIMAL;		}
	segsites:	 { return SEGSITES_TOKEN;	}
        prob:            { return PROB;                 }
	positions:	 { return POSITIONS_TOKEN;	}
	^\n		 { BEGIN(EAT); return BLANK_LINE;  }
	[\ \n\t]	 { /* EAT SPACES */	}
	.		 { return SCAN_ERROR;		}
}
.			 {  return SCAN_ERROR;		}
%%

yyscan_t msparser_setfile(FILE *f) {
  yyscan_t ms_scanner;
  struct yyguts_t *yyg;

  yylex_init(&ms_scanner);
  yyset_in(f, ms_scanner);

  yyg = ms_scanner;
  BEGIN(EAT);
  return ms_scanner;
}

#ifdef __cplusplus
 }
#endif
