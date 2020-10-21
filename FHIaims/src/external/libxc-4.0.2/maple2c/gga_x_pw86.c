/* 
  This file was generated automatically with /nfs/data-012/marques/software/source/libxc/svn/scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ../maple/gga_x_pw86.mpl
  Type of functional: work_gga_x
*/

void xc_gga_x_pw86_enhance
  (const xc_func_type *p,  xc_gga_work_x_t *r)
{
  double t1, t2, t3, t4, t5, t6, t7, t11;
  double t12, t14, t15, t19, t21, t25, t26, t27;
  double t29, t30, t31, t35, t42, t45, t46, t56;
  double t59;

  gga_x_pw86_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pw86_params * )(p->params);

  t1 = M_CBRT6;
  t2 = params->aa * t1;
  t3 = 0.31415926535897932385e1 * 0.31415926535897932385e1;
  t4 = cbrt(t3);
  t5 = t4 * t4;
  t6 = 0.1e1 / t5;
  t7 = r->x * r->x;
  t11 = t1 * t1;
  t12 = params->bb * t11;
  t14 = 0.1e1 / t4 / t3;
  t15 = t7 * t7;
  t19 = t3 * t3;
  t21 = params->cc / t19;
  t25 = 0.1e1 + t2 * t6 * t7 / 0.24e2 + t12 * t14 * t15 / 0.576e3 + t21 * t15 * t7 / 0.2304e4;
  r->f = pow(t25, 0.1e1 / 0.15e2);

  if(r->order < 1) return;

  t26 = r->f * r->f;
  t27 = t26 * t26;
  t29 = t27 * t27;
  t30 = t29 * t27 * t26;
  t31 = 0.1e1 / t30;
  t35 = t7 * r->x;
  t42 = t2 * t6 * r->x / 0.12e2 + t12 * t14 * t35 / 0.144e3 + t21 * t15 * r->x / 0.384e3;
  r->dfdx = t31 * t42 / 0.15e2;

  if(r->order < 2) return;

  t45 = 0.1e1 / t30 / t25;
  t46 = t42 * t42;
  t56 = t2 * t6 / 0.12e2 + t12 * t14 * t7 / 0.48e2 + 0.5e1 / 0.384e3 * t21 * t15;
  r->d2fdx2 = -0.14e2 / 0.225e3 * t45 * t46 + t31 * t56 / 0.15e2;

  if(r->order < 3) return;

  t59 = t25 * t25;
  r->d3fdx3 = 0.406e3 / 0.3375e4 / t30 / t59 * t46 * t42 - 0.14e2 / 0.75e2 * t45 * t42 * t56 + t31 * (t12 * t14 * r->x / 0.24e2 + 0.5e1 / 0.96e2 * t21 * t35) / 0.15e2;

  if(r->order < 4) return;


}

#define maple2c_order 3
#define maple2c_func  xc_gga_x_pw86_enhance
