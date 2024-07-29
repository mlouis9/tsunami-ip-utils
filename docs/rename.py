import os

def rename_and_remove_files(file_path):
    os.rename(f"{file_path}/tsunami_ip_utils.rst", f"{file_path}/index.rst")
    os.remove(f"{file_path}/modules.rst")

# Rename the tsunami_ip_utils.rst file to index.rst, and remove the modules.rst files
rename_and_remove_files("source/public_api")
rename_and_remove_files("source/private_api")

# Now move paths module to the root directory and rename it to paths.rst
os.rename("source/public_api/paths.rst", "source/paths.rst")
os.remove("source/private_api/paths.rst")