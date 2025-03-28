all	: urqmd.e

urqmd.e:
	$(MAKE) -C src
	cp -f src/urqmd.e ./
	$(MAKE) -C src distclean

clean:
	$(MAKE) -C src clean

distclean:
	$(MAKE) -C src distclean
	rm -f urqmd.e
