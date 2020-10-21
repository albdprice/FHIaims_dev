#
# Makefile for ATM
#
# Please note that ATM now depends on the GridXC and xmlf90 libraries
# ---------------------------------------------------------------
#
.SUFFIXES:
.SUFFIXES: .F90 .f .F .f90 .o .a 
#
default: atm
#
ARCH_MAKE?=arch.make
include $(ARCH_MAKE)
#
#
AUX_OBJS= auxf95.o
#
#--------------------------------------------------------------
dep:
	sfmakedepend --depend=obj  --modext=o *.f *.f90 *.F *.F90
#----------------------------------------------$$$$$
#
WXML=$(XMLF90_LIBS)
XC=$(GRIDXC_LIBS)
LIBS=$(WXML) $(XC)
#
MAIN_OBJS= atm.o  nucl_z.o   difnrl.o   difrel.o   dmixp.o   dsolv1.o   \
        pseudo.o  wf.o  hsc.o  excorr.o \
        tm2.o  ker.o  wrapup.o   splift.o   spliq.o  \
        velect.o   vionic.o   wtrans.o   v0pp.o  chg_mism.o  \
        dsolv2.o  etotal.o   ext.o  input.o  header.o orban.o   \
        potran.o   potrv.o potrvs.o  potrw.o   prdiff.o  string.o  root.o  \
        genrhs.o  logder.o ode.o denplot.o scilib.o  recipes.o \
	compat_params.o  pcc_exp.o  change_valence.o pswf_store.o \
	coreq.o get_unit.o gauleg.o divdif.o symbols.o m_uuid.o\
        prversion.o  write_psml.o pxf.o force_underflow.o dsolv2_single.o dsolv2_save_wf.o m_libxc_sxc_translation.o m_libxc_list.o m_siestaxc_list.o \
        handlers.o
#
OBJS= $(MAIN_OBJS) 
#
write_psml.o: INCFLAGS += $(XMLF90_INCFLAGS)
velect.o: INCFLAGS += $(GRIDXC_INCFLAGS)

ifeq ($(GRIDXC_USES_LIBXC),1)
velect.o: FPPFLAGS += -DLIBXC
endif

atm:	$(OBJS) $(AUX_OBJS) 
	$(FC) $(LDFLAGS) -o atm $(OBJS) $(AUX_OBJS) $(LIBS)
#
clean:	
	rm -f *.o *.mod atm *.a
#
# DO NOT DELETE THIS LINE - used by make depend
atm.o: charge.h param.h elecpot.h radial.h energy.h orbital.h
change_valence.o: ion.h param.h radial.h charge.h
chg_mism.o: radial.h nonlinear.h coeffs.h linear.h
compat_params.o: compat.h
coreq.o: radial.h charge.h
denplot.o: param.h radial.h charge.h
difnrl.o: elecpot.h ion.h param.h radial.h
difrel.o: param.h elecpot.h ion.h radial.h
dsolv1.o: charge.h radial.h energy.h orbital.h elecpot.h ion.h
dsolv2.o: ion.h elecpot.h param.h orbital.h energy.h radial.h charge.h
dsolv2_save_wf.o: charge.h radial.h energy.h orbital.h param.h elecpot.h ion.h
dsolv2_single.o: charge.h param.h elecpot.h ion.h orbital.h radial.h energy.h
etotal.o: orbital.h energy.h param.h
excorr.o: radial.h param.h
genrhs.o: nonlinear.h
header.o: compat.h version.h radial.h orbital.h param.h input.h
hsc.o: elecpot.h param.h ion.h energy.h radial.h orbital.h
input.o: param.h input.h version.h radial.h orbital.h compat.h
ker.o: elecpot.h ion.h param.h orbital.h energy.h radial.h
logder.o: energy.h ode_path.h ion.h radial.h orbital.h elecpot.h param.h
logder.o: ode_blk.h
m_libxc_sxc_translation.o: m_libxc_list.o m_siestaxc_list.o
nv0pp.o: coeffs.h
ode.o: ode_blk.h radial.h
orban.o: elecpot.h param.h ion.h radial.h
potran.o: plot.h
prversion.o: version.h
pseudo.o: elecpot.h param.h radial.h orbital.h
pswf_store.o: pseudowave.h radial.h ion.h
recipes.o: ode_path.h
symbols.o: set2.h
tm2.o: tm2_blk.h elecpot.h param.h ion.h orbital.h radial.h energy.h
v0pp.o: tm2_blk.h radial.h
velect.o: orbital.h energy.h radial.h param.h elecpot.h compat.h charge.h
velect.o: m_libxc_sxc_translation.o
vionic.o: charge.h orbital.h radial.h param.h ion.h
wf.o: energy.h radial.h orbital.h elecpot.h param.h ion.h
wrapup.o: charge.h elecpot.h param.h orbital.h radial.h version.h compat.h
wrapup.o: corecorr.h ion.h energy.h pseudowave.h
m_common.o: radial.h pseudowave.h orbital.h ion.h param.h corecorr.h charge.h
write_psml.o: m_common.o
write_psml.o: m_libxc_sxc_translation.o m_uuid.o
wtrans.o: plot.h
mtprng.o: m_uuid.o
