clean:
<<<<<<< HEAD
	@find ./ssat -name '__pycache__' -exec rm -rf {} \;
=======
	@find . -name '__pycache__' -exec rm -r {} \;
>>>>>>> 949afb8d030501ec029d5634fa4399a68a667bf2

.PHONY: clean
