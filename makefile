all:
	make -C csa_1.0/
	make -C radixSA/
	make -C kmismatch/

clean:
	make -C csa_1.0/ clean
	make -C radixSA/ clean
	make -C kmismatch/ clean
