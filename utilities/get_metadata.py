
import requests
import logging
from .classes.metadata import ToolMetadata

def get_metadata_from_biotools(biotoolsID):
    """
    Get metadata
    :param biotoolsID(str): biotoolsID from bio.tools
    :return:
      meta_dict(dict): dictionary of metadata from bio.tools.
    """
    params = {'format': 'json'}
    attrs = {'biotoolsID': f'"{biotoolsID}"'}
    r = requests.get(f"https://bio.tools/api/t/", params={**params, **attrs})
    biotools_dict = r.json()
    if biotools_dict['count'] != 1:
        logging.ERROR(f"bio.tools returned {biotools_dict['count']} results. Expected 1")
    return biotools_dict

def _handle_publication(publication_list):
    pub_list = []
    if publication_list:
        for publication in publication_list:
            pub_list.append({'identifier': publication.get('doi'), 'headline': publication['metadata'].get('title')})
        return pub_list
    else:
        return

def _handle_keywords(topics, functions):
    kew_words = [topic.get('uri') for topic in topics if topic.get('uri')]
    for function in functions:
        kew_words.extend([operation.get('uri') for operation in function['operation'] if operation.get('uri')])
    return kew_words

def _handle_credit(credit):
    creators = []
    contacts = []
    for entity in credit:
        if 'Developer' in entity['typeRole']:
            creators.append({'name': entity.get('name'), 'email': entity.get('email'), 'identifier': entity.get('orcidid')})
        if 'Primary contact' in entity['typeRole']:
            contacts.append({'name': entity.get('name'), 'email': entity.get('email'), 'identifier': entity.get('orcidid')})
    return {'creator': creators, 'contactPoint': contacts}


def pop_websites_and_repo(homepage, link, documentation):
    websites = []
    code_repo = None
    if homepage:
        websites.append({'name': None, 'description': 'homepage', 'URL': homepage})
    for entity in link:
        if entity['type'] == 'Repository':
            code_repo = {'name': None, 'URL': entity.get('url')}
        else:
            websites.append({'name': None, 'description': entity['type'], 'URL': entity.get('url')})
    for entity in documentation:
        websites.append({'name': None, 'description': 'documentation', 'URL': entity.get('url')})

    return {'WebSite': websites, 'codeRepository': code_repo}


def pop_metadata_template_from_biotools(biotoolsID):
    meta_dict = get_metadata_from_biotools(biotoolsID)
    meta_data = meta_dict['list'][0]
    tool_metadata = ToolMetadata(name=meta_data['name'],
                                 description=meta_data['description'],
                                 license=meta_data['license'],
                                 **pop_websites_and_repo(meta_data['homepage'], meta_data['link'], meta_data['documentation']),
                                 **_handle_credit(meta_data['credit']),
                                 publication=_handle_publication(meta_data.get('publication')),
                                 keywords=_handle_keywords(meta_data['topic'], meta_data['function']),
                                 # programmingLanguage=[],
                                 # datePublished='',
                                 # downloadURL=''
                                 )
    tool_metadata.extra.update({'biotools_id': biotoolsID})
    filename = f"{tool_metadata.name}-metadata.yaml"
    tool_metadata.mk_file(filename)
    return filename
