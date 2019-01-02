clean:
	@find ./ssat -name '__pycache__' -prune -exec rm -rf {} +

.PHONY: clean
