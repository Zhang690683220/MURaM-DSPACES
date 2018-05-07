default:
	make -C lib/CFD
	make -C src
clean: 
	make clean -C lib/CFD
	make clean -C src
        
distclean: 
	make distclean -C lib/CFD
	make distclean -C src
	rm -rf lib/{lib,include,share}
        

