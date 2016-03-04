
/*
 * Code to accompany the PHYS 462/562 Computational Biophysics course.
 * 
 * Copyright (c) 2009-2010 Luis R Cruz Cruz, Philadelphia, PA
 *
 * Please address questions to:
 *  Luis R Cruz Cruz
 *  Physics Department
 *  Drexel University
 *  3141 Chestnut St.
 *  Philadelphia, PA 19104
 */

#ifndef DEFS_H
#define DEFS_H

// DEFINES THE DIMENSIONALITY OF THE SYSTEM
#define NDIM 2
//#define NDIM 3

// DO CONFINED SYSTEM
#define CONFINED 0
//#define CONFINED 1

// DRAW VERLET NEIGHBORHOOD AROUND ONE PARTICLE
#define DRAWVERLETNEIGH 0
//#define DRAWVERLETNEIGH 1

// PUT IN TIME MEASURES TO QUANTIFY PERFORMANCE
#define BENCH_IT 0
//#define BENCH_IT 1

#define OUTPUT_FILES \
  ofstream params("results/params.txt"); \
  ofstream warn("results/warnings.txt"); \
  ofstream out_f("results/thermo.dat"); \
  ofstream out_p("results/totlinmom.dat"); \
  /* - CLEAN FILES TO BE APPENDED */ \
  ofstream out_gr_clean("results/gr.dat"); \
  ofstream out_r2t_clean("results/r2t.dat"); \
  ofstream out_vacf_clean("results/vacf.dat"); \
  ofstream out_tag_clean("results/tag.dat"); \
  ofstream out_rg("results/rg.dat"); \
  ofstream out_rg2("results/rg2.dat"); \
  ofstream out_eoe("results/eoe.dat"); \
  ofstream out_eoe_vec("results/eoe_vec.dat"); \
  ofstream out_eoe2("results/eoe2.dat"); \
  out_gr_clean << ""; \
  out_r2t_clean << ""; \
  out_vacf_clean << ""; \
  out_rg << ""; \
  out_rg2 << ""; \
  out_eoe << ""; \
  out_eoe_vec << ""; \
  out_eoe2 << ""; \
  out_gr_clean.close(); \
  out_r2t_clean.close(); \
  out_vacf_clean.close(); \
  out_tag_clean.close();  
 
#define CLEAN_OUTPUT_FILES \
  /* - CLEAN ALL FILES */ \
  ofstream params("results/params.txt"); \
  ofstream warn("results/warnings.txt"); \
  ofstream out_f("results/thermo.dat"); \
  ofstream out_p("results/totlinmom.dat"); \
  ofstream out_gr_clean("results/gr.dat"); \
  ofstream out_r2t_clean("results/r2t.dat"); \
  ofstream out_vacf_clean("results/vacf.dat"); \
  ofstream out_tag_clean("results/tag.dat"); \
  ofstream out_rg("results/rg.dat"); \
  ofstream out_rg2("results/rg2.dat"); \
  ofstream out_eoe("results/eoe.dat"); \
  ofstream out_eoe_vec("results/eoe_vec.dat"); \
  ofstream out_eoe2("results/eoe2.dat"); \
  params << ""; \
  warn << ""; \
  out_f << ""; \
  out_p << ""; \
  out_rg << ""; \
  out_rg2 << ""; \
  out_eoe << ""; \
  out_eoe_vec << ""; \
  out_eoe2 << ""; \
  out_gr_clean << ""; \
  out_r2t_clean << ""; \
  out_vacf_clean << ""; \
  params.close(); \
  warn.close(); \
  out_f.close(); \
  out_p.close(); \
  out_gr_clean.close(); \
  out_r2t_clean.close(); \
  out_vacf_clean.close(); \
  out_tag_clean.close();  \
  out_rg.close();  \
  out_rg2.close();  \
  out_eoe.close();  \
  out_eoe_vec.close();  \
  out_eoe2.close(); 


#define SET_SEED \
  long seed; \
  long time_sub = 1234234670; \
  /* USE TIME AS SEED */ \
  seed = - long( time(NULL) - time_sub ); \
  /* seed = -3;*/ \
  /* SEED THE GENERATOR */ \
  srand(seed);


#endif // DEFS_H

