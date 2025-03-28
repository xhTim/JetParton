all	: osc2u.e

osc2u.e:
	$(MAKE) -C src
	cp -f src/osc2u.e ./

clean:
	$(MAKE) -C src clean

distclean:
	$(MAKE) -C src distclean
	rm -f osc2u.e
