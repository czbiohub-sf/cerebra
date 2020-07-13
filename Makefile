help:
	@echo "lint - check code style with flake8"
	@echo "test - run tests only"
	@echo "coverage - run tests and check code coverage"
	@echo "install - Install requirements "

test:
	py.test

coverage:
	coverage run --source cerebra --omit="*/test*" --module py.test
	coverage report --show-missing

lint:
	flake8 --exclude docs cerebra --exit-zero

install:
	pip install -r requirements.txt
	pip install .
