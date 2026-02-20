import pandas as pd

def clean_text(text):
    '''
    Remove nonstandard characters except mutation arrows.
    '''
    if pd.isna(text) or not isinstance(text, str):
        return text
    # Replace non-breaking hyphen with regular hyphen
    text = text.replace('‑', '-')

    chars = []
    for character in text:
        # Remove nonstandard characters except for mutation arrows
        if ord(character) < 128 or character in ['→', '←']:
            chars.append(character)

    return ''.join(chars)