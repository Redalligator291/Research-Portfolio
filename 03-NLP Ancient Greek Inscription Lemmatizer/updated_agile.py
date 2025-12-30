"""
Original agile.py by Jasper bos modified to solve problem with Stanza tokenizer.
Modified by Silvia Stopponi on 15/05/2023 (lines 39ss)
"""
from cltk.alphabet.grc import grc  # For normalizing texts
import re
import stanza
import pickle
from Levenshtein import distance
from stanza.models.common.doc import Document
import unicodedata as ud
from unidecode import unidecode
​
def remove_accents(accenty):
    nfkd_form = ud.normalize('NFKD', accenty)
    return "".join([c for c in nfkd_form if not ud.combining(c)])
​
def is_capital(wordy):
    if wordy[0].isupper():
        return True
    return False
​
def capitalization(wordy):
    return wordy[0].upper() + wordy[1:]
​
​
def distance_accents(entry, predicted):
    dist = distance(entry, predicted)
    max_index = min(len(entry), len(predicted))
​
    for index in range(max_index):
        if entry[index] != predicted[index] and remove_accents(entry[index]) == remove_accents(predicted[index]):
            dist -= 1
    return dist
        
​
def lemmatize(text, use_lexicon=True):
    """
    Lemmatizes (and tokenizes) Ancient Greek inscriptions given custom rules,
    a custom trained model in Stanza and a lexicon lookup.
​
    :param text: inscription text to be lemmatized
    :param use_lexicon: bool to enable or disable the correction by the lexicon lookup
    :return: Stanza Document containing id, text, lemma, start_char and end_char annotations
    """
​
    # Handle extra-alphabetical characters
    original_text = grc.normalize_grc(text)
    original_text = re.sub('[|∣·∙∶:,.⁝⋮⁞⁙“”]+', ' \g<0> ', original_text)  # Pre-tokenize
    original_text = re.sub(' +', ' ', original_text)
    original_text = re.sub('\n+', '\n', original_text)  # Sentence tokenization not supported
​
    # Add custom rules (chars have been normalized)
    processed_text = re.sub('(?<!\s)[Ϝϝh](?!\s)', '', original_text)  # [Ϝϝh] within token
    processed_text = re.sub('(?<=\s)[Ϝϝh](?!\s)(?=.)', '', processed_text)  # [Ϝϝh] begin of token
    processed_text = re.sub('(?<=.)(?<!\s)[Ϝϝh](?=\s)', '', processed_text)  # [Ϝϝh] end of token
    processed_text = re.sub('(κς)|(κσ)|(χς)|(χσ)', 'ξ', processed_text)
    processed_text = re.sub('(Κς)|(Κσ)|(Χσ)|(Χς)', 'Ξ', processed_text)
    processed_text = re.sub('(φς)|(φσ)', 'ψ', processed_text)
    processed_text = re.sub('(Φς)|(Φσ)', 'Ψ', processed_text)
    processed_text = re.sub(' [|∣·∙∶:,.⁝⋮⁞⁙“”]+', '', processed_text)
​
    # old code:
    # lemma_nlp = stanza.Pipeline(lang='grc', processors='tokenize,lemma', tokenize_no_ssplit=True,
    #                             lemma_model_path='grc_agile_lemmatizer.pt', verbose=False)
    # token_nlp = stanza.Pipeline(lang='grc', processors='tokenize', tokenize_no_ssplit=True, verbose=False)
    # new  code:
    lemma_nlp = stanza.Pipeline(lang='grc', processors='tokenize,lemma', tokenize_pretokenized=True,
                                lemma_model_path='grc_agile_lemmatizer.pt', verbose=False)
    token_nlp = stanza.Pipeline(lang='grc', processors='tokenize', tokenize_pretokenized=True, verbose=False)
    token_dict = token_nlp(original_text).to_dict()[0]  # Dict for all tokens (lemmas to be inserted)
    lemma_dict = lemma_nlp(processed_text).to_dict()[0]  # Dict for lemmas given by model
    lexicon = pickle.load(open("lexicon.p", "rb"))
​
​
    token_dict_original = token_nlp(original_text).to_dict()[0]  # Dict for all tokens (lemmas to be inserted)
​
​
    # Add lemmas to token dict
    lemma_i = 0
​
    for token_i, token in enumerate(token_dict):
        if re.search('[|∣·∙∶:,.⁝⋮⁞⁙“”]+', token['text']):  # Custom lemmatization
            token_dict[token_i]['lemma'] = token['text']
        else:
            try:
                predicted = lemma_dict[lemma_i]['lemma']  # Lemmatization by model
            except KeyError:  # No lemma
                predicted = ""
            # Handle lexicon correction
            if use_lexicon:
​
                if predicted != "" and predicted not in lexicon:
                    
                    lowest = distance_accents(lexicon[0], predicted)
                    closest = lexicon[0]
​
                    for entry in lexicon[1:]:
                        dist = distance_accents(entry, predicted)
                        if dist == 0:  # Speed optimisation
                            closest = entry
                            lowest = 0
                            break
                        if dist < lowest:
                            lowest = dist
                            closest = entry
                            print("lowest distance:", lowest, "closest word:", closest, "predicted word:", predicted)
                        
​
                    token_dict[token_i]['lemma'] = closest
                else:
                    token_dict[token_i]['lemma'] = predicted
            else:
                token_dict[token_i]['lemma'] = predicted
            
​
            
            if is_capital(token_dict_original[token_i]['text']):
                token_dict[token_i]['lemma'] = capitalization(token_dict[token_i]['lemma'])
            lemma_i += 1
​
    return Document([token_dict])
​
#through lexicon, minimum distance, find all with minimum distance in a new list. If list size = 1, continue. 
# Else, needs remove_accent for each item in list.
#           ex if list = [jello, hEllo], list will be converted without accents. then calculate distance again.
#           Also, remove accents on actual output word
#           Do edit distance again between each word of the modified list with the output word of the model
#           output the least edit distance