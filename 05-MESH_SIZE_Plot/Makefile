.PHONY: all

all:
	git clone --filter=blob:none --sparse https://github.com/chapuisl/CFD_TOOLS.git temp_repo
	cd temp_repo && git sparse-checkout set 05-MESH_SIZE_Plot
	mv temp_repo/05-MESH_SIZE_Plot .
	rm -rf temp_repo
	rm -f Makefile
