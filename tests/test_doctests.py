import doctest
import os
import importlib
import pkgutil

def find_and_run_doctests(package_name):
    package = importlib.import_module(package_name)
    print(f"Testing package: {package_name}")
    print('-'*len(f"Testing package: {package_name}"))
    errors = 0
    failures = 0

    for importer, modname, ispkg in pkgutil.walk_packages(package.__path__, package.__name__ + '.'):
        try:
            module = importlib.import_module(modname)
            print(f"Testing module: {modname}")
            result = doctest.testmod(module)
            failures += result.failed
            if result.failed:
                print(f"Failed: {result.failed}, Attempted: {result.attempted}")
        except Exception as e:
            print(f"Failed to test module {modname}: {str(e)}")

    print(f"\nFinished testing. Failures: {failures}")

# Replace 'mypackage' with your actual package name.
if __name__ == "__main__":
    find_and_run_doctests('tsunami_ip_utils')