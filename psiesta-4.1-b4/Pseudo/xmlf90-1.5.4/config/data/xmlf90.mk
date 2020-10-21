# XMLF90 helper mk
#
# In order to hide the structure of this hierarchy, this file
# will export the needed symbols for a client.
#
# All you need is to define XMLF90_ROOT in your makefile
# to point to your installation of the XMLF90 library.
#
#----------------------------------------------------------
#
ifdef XMLF90_H__

$(info multiple inclusions of xmlf90.mk ...)

else

XMLF90_H__ = 1

XMLF90_ROOT_BUILD = /home/lleblanc/src/psiesta-4.1-b4/Pseudo/xmlf90-1.5.4
XMLF90_INCFLAGS   = -I$(XMLF90_ROOT_BUILD)/include 
XMLF90_LIBS       = $(XMLF90_ROOT_BUILD)/lib/libxmlf90.a

endif
