from docutils import nodes
from docutils.parsers.rst import roles
from tsunami_ip_utils import config

def global_param_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    # Fetch the global parameter value using dictionary access
    value = getattr(config, text, 'undefined')
    node = nodes.literal('', str(value), **options)  # Create a node with the value
    return [node], []

def setup(app):
    roles.register_local_role('globalparam', global_param_role)