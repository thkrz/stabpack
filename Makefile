STABPACK = ssat/ssatlib/stabpack

all: stabpack

stabpack:
	@echo BUILD STABPACK
	@$(MAKE) --no-print-directory -C ${STABPACK}

clean:
	@echo cleaning...
	@find ./ssat -name '__pycache__' -prune -exec rm -rf {} +
	@$(MAKE) --no-print-directory -C ${STABPACK} clean

.PHONY: clean
