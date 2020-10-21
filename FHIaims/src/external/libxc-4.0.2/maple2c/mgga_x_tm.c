/* 
  This file was generated automatically with /nfs/data-012/marques/software/source/libxc/svn/scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ../maple/mgga_x_tm.mpl
  Type of functional: work_mgga_x
*/

static void 
xc_mgga_x_tm_enhance(const xc_func_type *pt, xc_mgga_work_x_t *r)
{
  double t1, t2, t3, t4, t5, t6, t7, t10;
  double t12, t13, t15, t17, t18, t19, t20, t21;
  double t22, t23, t24, t25, t27, t29, t30, t33;
  double t34, t38, t45, t46, t47, t50, t52, t55;
  double t65, t66, t68, t69, t73, t77, t78, t81;
  double t85, t87, t94, t95, t98, t100, t101, t102;
  double t103, t104, t105, t107, t110, t116, t121, t122;
  double t123, t125, t126, t129, t130, t131, t132, t135;
  double t140, t142, t144, t153, t156, t159, t163, t165;
  double t184, t185, t199, t200, t203, t205, t210, t212;
  double t213, t215, t216, t220, t221, t240, t243;


  t1 = M_CBRT6;
  t2 = 0.31415926535897932385e1 * 0.31415926535897932385e1;
  t3 = POW_1_3(t2);
  t4 = t3 * t2;
  t5 = t1 * t4;
  t6 = r->t * r->t;
  t7 = 0.1e1 / t6;
  t10 = t2 * t2;
  t12 = 0.1e1 / t6 / r->t;
  t13 = t10 * t12;
  t15 = 0.27e2 / 0.50e2 * t5 * t7 + 0.729e3 / 0.250e3 * t13;
  t17 = 0.1e1 + 0.243e3 / 0.250e3 * t13;
  t18 = t17 * t17;
  t19 = 0.1e1 / t18;
  t20 = t15 * t19;
  t21 = t3 * t3;
  t22 = 0.1e1 / t21;
  t23 = t1 * t22;
  t24 = r->x * r->x;
  t25 = t23 * t24;
  t27 = t1 * t1;
  t29 = t27 / t4;
  t30 = t24 * t24;
  t33 = 0.1e1 + 0.15045488888888888889e0 * t25 + 0.26899490462262948000e-2 * t29 * t30;
  t34 = pow(t33, 0.1e1 / 0.5e1);
  t38 = t27 * t21;
  t45 = 0.1e1 + 0.63943327777777777778e-1 * t25 - 0.5e1 / 0.9e1 * (0.14554132000000000000e0 * r->t + 0.25633760400000000000e0 * t38 + 0.11867481666666666667e-1 * t24) * t1 * t22;
  t46 = t34 * t34;
  t47 = 0.1e1 / t46;
  t50 = 0.1e1 / t34 + 0.7e1 / 0.9e1 * t45 * t47;
  t52 = 0.1e1 - t20;
  t55 = (0.10e2 / 0.81e2 + 0.25e2 / 0.8748e4 * t25) * t1;
  t65 = (r->t - t24 / 0.8e1) * t1 * t22 / 0.4e1 - 0.9e1 / 0.20e2 + t25 / 0.36e2;
  t66 = t65 * t65;
  t68 = t65 * t27;
  t69 = 0.1e1 / r->t;
  t73 = 0.1e1 - 0.3e1 / 0.10e2 * t38 * t69;
  t77 = 0.1e1 + 0.5e1 / 0.12e2 * t55 * t22 * t24 + 0.292e3 / 0.405e3 * t66 - 0.73e2 / 0.225e3 * t68 * t21 * t69 * t73;
  t78 = pow(t77, 0.1e1 / 0.10e2);
  r->f = t20 * t50 + t52 * t78;

  if(r->order < 1) return;

  r->dfdrs = 0.0e0;
  t81 = 0.1e1 / t34 / t33;
  t85 = t29 * t24 * r->x;
  t87 = 0.30090977777777777778e0 * t23 * r->x + 0.10759796184905179200e-1 * t85;
  t94 = 0.1e1 / t46 / t33;
  t95 = t45 * t94;
  t98 = -t81 * t87 / 0.5e1 + 0.89211550411522633749e-1 * t23 * r->x * t47 - 0.14e2 / 0.45e2 * t95 * t87;
  t100 = t78 * t78;
  t101 = t100 * t100;
  t102 = t101 * t101;
  t103 = t102 * t78;
  t104 = 0.1e1 / t103;
  t105 = t52 * t104;
  t107 = t22 * r->x;
  t110 = t65 * t1;
  t116 = 0.125e3 / 0.52488e5 * t85 + 0.5e1 / 0.6e1 * t55 * t107 - 0.73e2 / 0.7290e4 * t110 * t107 + 0.73e2 / 0.5400e4 * r->x * t69 * t73;
  r->dfdx = t20 * t98 + t105 * t116 / 0.10e2;
  t121 = t6 * t6;
  t122 = 0.1e1 / t121;
  t123 = t10 * t122;
  t125 = -0.27e2 / 0.25e2 * t5 * t12 - 0.2187e4 / 0.250e3 * t123;
  t126 = t125 * t19;
  t129 = 0.1e1 / t18 / t17;
  t130 = t15 * t129;
  t131 = t50 * t10;
  t132 = t131 * t122;
  t135 = t23 * t47;
  t140 = -t126 - 0.729e3 / 0.125e3 * t130 * t123;
  t142 = t110 * t22;
  t144 = t69 * t73;
  t153 = 0.146e3 / 0.405e3 * t142 - 0.73e2 / 0.150e3 * t144 + 0.73e2 / 0.225e3 * t68 * t21 * t7 * t73 - 0.73e2 / 0.125e3 * t110 * t4 * t12;
  r->dfdt = t126 * t50 + 0.729e3 / 0.125e3 * t130 * t132 - 0.62888224691358024691e-1 * t20 * t135 + t140 * t78 + t105 * t153 / 0.10e2;
  r->dfdu = 0.0e0;

  if(r->order < 2) return;

  r->d2fdrs2 = 0.0e0;
  t156 = t33 * t33;
  t159 = t87 * t87;
  t163 = t29 * t24;
  t165 = 0.30090977777777777778e0 * t23 + 0.32279388554715537600e-1 * t163;
  t184 = t52 / t103 / t77;
  t185 = t116 * t116;
  r->d2fdx2 = t20 * (0.6e1 / 0.25e2 / t34 / t156 * t159 - t81 * t165 / 0.5e1 + 0.89211550411522633749e-1 * t135 - 0.71369240329218107000e-1 * t23 * r->x * t94 * t87 + 0.98e2 / 0.225e3 * t45 / t46 / t156 * t159 - 0.14e2 / 0.45e2 * t95 * t165) - 0.9e1 / 0.100e3 * t184 * t185 + t105 * (0.1397e4 / 0.116640e6 * t163 + 0.5e1 / 0.6e1 * t55 * t22 - 0.73e2 / 0.7290e4 * t142 + 0.73e2 / 0.5400e4 * t144) / 0.10e2;
  t199 = 0.1e1 / t121 / r->t;
  t200 = t10 * t199;
  t203 = (0.81e2 / 0.25e2 * t5 * t122 + 0.4374e4 / 0.125e3 * t200) * t19;
  t205 = t125 * t129;
  t210 = t18 * t18;
  t212 = t15 / t210;
  t213 = t10 * t10;
  t215 = t121 * t121;
  t216 = 0.1e1 / t215;
  t220 = t130 * t1;
  t221 = t22 * t47;
  t240 = t140 * t104;
  t243 = t153 * t153;
  r->d2fdt2 = t203 * t50 + 0.1458e4 / 0.125e3 * t205 * t132 - 0.12577644938271604938e0 * t126 * t135 + 0.1594323e7 / 0.31250e5 * t212 * t50 * t213 * t216 - 0.35726160176503976590e2 * t220 * t221 * t122 - 0.2916e4 / 0.125e3 * t130 * t131 * t199 - 0.36676412640000000000e0 * t220 * t221 * t123 + (-t203 - 0.1458e4 / 0.125e3 * t205 * t123 - 0.1594323e7 / 0.31250e5 * t212 * t213 * t216 + 0.2916e4 / 0.125e3 * t130 * t200) * t78 + t240 * t153 / 0.5e1 - 0.9e1 / 0.100e3 * t184 * t243 + t105 * (0.73e2 / 0.810e3 * t29 + 0.73e2 / 0.75e2 * t7 * t73 - 0.73e2 / 0.250e3 * t12 * t27 * t21 - 0.146e3 / 0.225e3 * t68 * t21 * t12 * t73 + 0.292e3 / 0.125e3 * t110 * t4 * t122) / 0.10e2;
  r->d2fdu2 = 0.0e0;
  r->d2fdrsx = 0.0e0;
  r->d2fdrst = 0.0e0;
  r->d2fdrsu = 0.0e0;
  r->d2fdxt = t126 * t98 + 0.729e3 / 0.125e3 * t130 * t98 * t10 * t122 + 0.25155289876543209877e-1 * t20 * t1 * t22 * t94 * t87 + t240 * t116 / 0.10e2 - 0.9e1 / 0.100e3 * t184 * t116 * t153 + t105 * (-0.73e2 / 0.29160e5 * t29 * r->x - 0.73e2 / 0.5400e4 * r->x * t7 * t73 + 0.73e2 / 0.18000e5 * r->x * t12 * t38) / 0.10e2;
  r->d2fdxu = 0.0e0;
  r->d2fdtu = 0.0e0;

  if(r->order < 3) return;


}

#define maple2c_order 3
#define maple2c_func  xc_mgga_x_tm_enhance
