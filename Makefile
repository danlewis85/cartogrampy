##
# Cluster App
#
# @file
# @version 0.1

.DEFAULT_GOAL := install


#-- General ---------------------------------------------------------------------
.PHONY: pipenv install

install:
	@printf "Setting up python environment...\n"
	@pip install --upgrade -r requirements.txt
	@printf "done..\n"

pipenv:
	@printf "Setting up pipenv environment...\n"
	@pip install --upgrade pipenv
	@pipenv install -r requirements.txt
	@printf "done..\n"


# end
