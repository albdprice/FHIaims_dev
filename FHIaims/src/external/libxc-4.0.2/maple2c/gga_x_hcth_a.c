/* 
  This file was generated automatically with /nfs/data-012/marques/software/source/libxc/svn/scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ../maple/gga_x_hcth_a.mpl
  Type of functional: work_gga_x
*/

void xc_gga_x_hcth_a_enhance
  (const xc_func_type *p,  xc_gga_work_x_t *r)
{
  double t1, t2, t3, t4, t6, t7, t8, t9;
  double t10, t13, t16, t17, t19, t23, t28, t29;
  double t30, t33, t37, t38, t40, t50, t55, t58;
  double t61, t62, t67, t77, t85, t90;


  t1 = M_CBRT3;
  t2 = t1 * t1;
  t3 = M_CBRT4;
  t4 = t2 * t3;
  t6 = cbrt(0.1e1 / 0.31415926535897932385e1);
  t7 = 0.1e1 / t6;
  t8 = r->x * r->x;
  t9 = t7 * t8;
  t10 = log(r->x + sqrt(r->x * r->x + 0.1e1));
  t13 = 0.1e1 + 0.252e-1 * r->x * t10;
  t16 = t13 * t13;
  t17 = 0.1e1 / t16;
  t19 = -0.251173e1 / t13 + 0.37198333333333333333e1 * t17;
  r->f = 0.109878e1 + 0.93333333333333333332e-3 * t4 * t9 * t19;

  if(r->order < 1) return;

  t23 = t7 * r->x;
  t28 = t8 + 0.1e1;
  t29 = sqrt(t28);
  t30 = 0.1e1 / t29;
  t33 = 0.252e-1 * t10 + 0.252e-1 * r->x * t30;
  t37 = 0.1e1 / t16 / t13;
  t38 = t37 * t33;
  t40 = 0.251173e1 * t17 * t33 - 0.74396666666666666666e1 * t38;
  r->dfdx = 0.18666666666666666666e-2 * t4 * t23 * t19 + 0.93333333333333333332e-3 * t4 * t9 * t40;

  if(r->order < 2) return;

  t50 = t33 * t33;
  t55 = 0.1e1 / t29 / t28;
  t58 = 0.504e-1 * t30 - 0.252e-1 * t8 * t55;
  t61 = t16 * t16;
  t62 = 0.1e1 / t61;
  t67 = -0.502346e1 * t37 * t50 + 0.251173e1 * t17 * t58 + 0.22319000000000000000e2 * t62 * t50 - 0.74396666666666666666e1 * t37 * t58;
  r->d2fdx2 = 0.18666666666666666666e-2 * t4 * t7 * t19 + 0.37333333333333333332e-2 * t4 * t23 * t40 + 0.93333333333333333332e-3 * t4 * t9 * t67;

  if(r->order < 3) return;

  t77 = t50 * t33;
  t85 = t28 * t28;
  t90 = -0.1008e0 * t55 * r->x + 0.756e-1 * t8 * r->x / t29 / t85;
  r->d3fdx3 = 0.55999999999999999998e-2 * t4 * t7 * t40 + 0.55999999999999999998e-2 * t4 * t23 * t67 + 0.93333333333333333332e-3 * t4 * t9 * (0.1507038e2 * t62 * t77 - 0.1507038e2 * t38 * t58 + 0.251173e1 * t17 * t90 - 0.89276000000000000000e2 / t61 / t13 * t77 + 0.66957000000000000000e2 * t62 * t33 * t58 - 0.74396666666666666666e1 * t37 * t90);

  if(r->order < 4) return;


}

#define maple2c_order 3
#define maple2c_func  xc_gga_x_hcth_a_enhance
