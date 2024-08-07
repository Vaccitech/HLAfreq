# noxfile.py
import nox


@nox.session
def lint(session):
    session.install("flake8")
    session.run("flake8", "src", "--count", "--select=E9,F63,F7,F82", "--show-source", "--statistics")
    session.run("flake8", "src", "--count", "--exit-zero", "--max-complexity=10", "--max-line-length=127", "--statistics")

@nox.session
def bandit(session):
    session.install("bandit")
    session.run("bandit", "-r", "src")

@nox.session(python=["3.10", "3.11", "3.12"])
def tests(session):
    session.install("pytest", "pytest-cov")
    session.install(".")
    session.run("pytest", "--cov=HLAfreq")
