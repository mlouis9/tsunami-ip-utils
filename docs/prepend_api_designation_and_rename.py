import os

def prepend_text(file_path, text):
    with open(file_path, 'r+') as file:
        content = file.read()
        # Check if the text is already at the start of the file
        if not content.startswith(text.rstrip('\r\n')):
            file.seek(0, 0)  # Rewind file to the beginning
            file.write(text.rstrip('\r\n') + '\n\n' + content)  # Prepend text

def update_rst_files(directory, label):
    for filename in os.listdir(directory):
        if filename.endswith(".rst"):
            file_path = os.path.join(directory, filename)
            prepend_text(file_path, f".. _{label}:")

def rename_and_remove_files(file_path):
    os.rename(f"{file_path}/tsunami_ip_utils.rst", f"{file_path}/index.rst")
    os.remove(f"{file_path}/modules.rst")

# Update public API files
update_rst_files("source/public_api", "public_api")

# Update private API files
update_rst_files("source/private_api", "private_api")

# Rename the tsunami_ip_utils.rst file to index.rst, and remove the modules.rst files
rename_and_remove_files("source/public_api")
rename_and_remove_files("source/private_api")

