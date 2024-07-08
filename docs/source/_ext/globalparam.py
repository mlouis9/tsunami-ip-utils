from docutils import nodes
from docutils.parsers.rst import roles
import tsunami_ip_utils  # import your module where the global parameter is defined

def global_param_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    # Fetch the global parameter value using dictionary access
    value = tsunami_ip_utils.config.get(text, 'undefined')  # Use dict.get() to handle missing keys gracefully
    node = nodes.literal('', str(value), **options)  # Create a node with the value
    return [node], []

def setup(app):
    roles.register_local_role('globalparam', global_param_role)