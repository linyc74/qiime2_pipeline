"""
LEfSe is terrible at treating non-alphanumeric characters
This is a mapping of characters to their replacements
"""
ORIGINAL_TO_NEW = [
    (' ', 'BLANK'),
    ('-', 'HYPHEN'),
    ('.', 'DOT'),
    ('[', 'OPENBRACKET'),
    (']', 'CLOSEBRACKET'),
    ('(', 'OPENPAREN'),
    (')', 'CLOSEPAREN'),
    ('{', 'OPENBRACE'),
    ('}', 'CLOSEBRACE'),
    ('/', 'SLASH'),
]
