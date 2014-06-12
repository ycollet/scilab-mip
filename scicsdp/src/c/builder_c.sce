// This file is released into the public domain

build_openmp = %f;

src_c_path = get_absolute_file_path("builder_c.sce");

// Use this line to specify options for the C compiler.  You'll probably
// want to turn on optimizations. You may also have to use some of the 
// following flags:
//
//  -DCAPSBLAS         if BLAS routine names are capitalized.
//  -DCAPSLAPACK       if LAPACK routine names are capitalized.
//  -DNOUNDERBLAS      if BLAS routine names have no underscore.
//  -DNOUNDERLAPACK    if LAPACK routine names have no underscore.
//  -DBIT64            For I32LP64 systems.
//  -DNOSHORTS         Allow for (LP) blocks of more than 65535 variables.
//  -DUSEOPENMP        Build an OpenMP parallel version.
//  -DSETNUMTHREADS    Work with OpenMP aware BLAS.  
//  -DUSESIGTERM       Terminate nicely at the end of the next iteration
//                     after receiving a SIGTERM signal
//  -DUSEGETTIME       Use ANSI C gettime() routine to determine clock
//                     time used in different parts of the code.
//  -DUSEATLAS         Turns on some special code for use with the ATLAS BLAS.

Symbols = ['sdp','triu','store_packed','store_unpacked','alloc_mat_packed',...
           'free_mat_packed','structnnz','actnnz','bandwidth',...
           'qreig','sort_entries','norm2','norm1','norminf',...
           'Fnorm','Knorm','mat1norm','matinfnorm','calc_pobj',...
           'calc_dobj','trace_prod','linesearch','pinfeas',...
           'dinfeas','dimacserr3','op_a','op_at','makefill',...
           'op_o','addscaledmat','zero_mat','add_mat','sym_mat',...
           'make_i','copy_mat','mat_mult','mat_multspa','mat_multspb',...
           'mat_multspc','mat_mult_raw','mat_mult_rawatlas','matvec',...
           'alloc_mat','free_mat','initparams','initsoln','trans',...
           'chol_inv','chol','solvesys','user_exit','read_sol',...
           'read_prob','write_prob','write_sol','free_prob',...
           'easy_sdp','sci_easy_sdp','tweakgap','bisect_'];

Files = ['add_mat.c','calc_dobj.c','copy_mat.c','freeprob.c','linesearch.c',...
         'mat_mult.c','norms.c','op_o.c','qreig.c','sdp.c','sym_mat.c',...
         'user_exit.c','zero_mat.c','addscaledmat.c','calc_pobj.c','easysdp.c','sci_easysdp.c',...
         'initparams.c','makefill.c','mat_multsp.c','op_a.c','packed.c',...
         'readprob.c','solvesys.c','trace_prod.c','writeprob.c','allocmat.c',...
         'chol.c','Fnorm.c','initsoln.c','make_i.c','matvec.c','op_at.c',...
         'psd_feas.c','readsol.c','sortentries.c','tweakgap.c','writesol.c',...
	 'overload_printf.c'];

if getos()~='Windows' then
  if build_openmp then
    LDFLAGS = '-lgomp -Wl,--wrap,printf';
  else
    LDFLAGS = '-lgomp -Wl,--wrap,printf';
  end
else
  LDFLAGS = '';
end

if build_openmp then
  CFLAGS = '-g -Wall -DNOSHORTS -DBIT64 -DUSEOPENMP -DSETNUMTHREADS -I' + src_c_path;
else
  CFLAGS = '-g -Wall -DNOSHORTS -DBIT64 -I' + src_c_path;
end

tbx_build_src(Symbols, Files, 'c', src_c_path, '', LDFLAGS, CFLAGS);

clear tbx_build_src;
clear src_c_path;
clear CFLAGS;
clear LDFLAGS;
