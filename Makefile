# To use this Makefile, get a copy of my SF Release Tools
# git clone git://git.code.sf.net/p/sfreleasetools/code sfreleasetools
# And point the environment variable RELEASETOOLS to the checkout
PNAME=antenna-optimizer
PYNAME=antenna_optimizer
ifeq (,${RELEASETOOLS})
    RELEASETOOLS=../releasetools
endif
LASTRELEASE:=$(shell $(RELEASETOOLS)/lastrelease -n)
PYF=antenna_model.py coaxmodel.py folded_3ele.py folded_bc.py \
    folded_bigrefl.py folded.py hb9cv.py logper.py statstool.py tl.py \
    transmission.py
VERSIONPY=$(PYNAME)/Version.py
VERSION=$(VERSIONPY)
README=README.rst
SRC=Makefile setup.py $(PYF:%.py=$(PYNAME)/%.py) \
    MANIFEST.in $(README) README.html

USERNAME=schlatterbeck
PROJECT=$(PNAME)
PACKAGE=$(PNAME)
CHANGES=changes
NOTES=notes

all: $(VERSION)

$(VERSION): $(SRC)

dist: all
	$(PYTHON) setup.py sdist --formats=gztar,zip

clean:
	rm -f $(PYNAME)/Version.py default.css            \
	      README.aux README.dvi README.log README.out \
	      README.tex
	rm -rf $(CLEAN)

include $(RELEASETOOLS)/Makefile-sf
