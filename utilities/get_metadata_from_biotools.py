
import requests
import logging
# from .classes.metadata import ToolMetadata

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
        logging.error(f"bio.tools returned {biotools_dict['count']} results. Expected 1")
    return biotools_dict

def _handle_publication(publication_list):
    pub_list = []
    if publication_list:
        for publication in publication_list:
            pub_dict = {'identifier': publication.get('doi')}
            if publication.get('metadata'):
                pub_dict['headline'] = publication['metadata'].get('title')
            pub_list.append(pub_dict)
        return pub_list
    else:
        return

def _handle_keywords(topics, functions):
    key_words = [topic.get('uri') for topic in topics if topic.get('uri')]
    for function in functions:
        key_words.extend([operation.get('uri') for operation in function['operation'] if operation.get('uri')])
    return key_words if key_words else None

def _handle_credit(credit):
    if credit:
        creators = []
        contacts = []
        for entity in credit:
            if 'Developer' in entity['typeRole']:
                creators.append({'name': entity.get('name'), 'email': entity.get('email'), 'identifier': entity.get('orcidid')})
            if 'Primary contact' in entity['typeRole']:
                contacts.append({'name': entity.get('name'), 'email': entity.get('email'), 'identifier': entity.get('orcidid')})
        return_dict = {'creator': creators if creators else None, 'contactPoint': contacts if contacts else None}
    else:
        return_dict = {'creator': None, 'contactPoint': None}
    return return_dict


def pop_websites_and_repo(homepage, link, documentation):
    websites = []
    code_repo = {'name': None, 'URL': None}
    if homepage:
        websites.append({'name': None, 'description': 'homepage', 'URL': homepage})
    for entity in link:
        if entity['type'] == 'Repository':
            code_repo = {'name': None, 'URL': entity.get('url')}
        else:
            websites.append({'name': None, 'description': entity['type'], 'URL': entity.get('url')})
    for entity in documentation:
        websites.append({'name': None, 'description': 'documentation', 'URL': entity.get('url')})

    return {'WebSite': websites if websites else None, 'codeRepository': code_repo}


def make_tool_metadata_kwargs_from_biotools(biotools_id):
    meta_dict = get_metadata_from_biotools(biotools_id)
    meta_data = meta_dict['list'][0]
    tool_kwargs = {}
    tool_kwargs['name'] = meta_data['name']
    tool_kwargs['description'] = meta_data['description']
    tool_kwargs['license'] = meta_data['license']
    tool_kwargs['publication'] = _handle_publication(meta_data.get('publication'))
    tool_kwargs['programmingLanguage'] = meta_data.get('language', None)
    tool_kwargs['keywords'] = _handle_keywords(meta_data['topic'], meta_data['function'])
    tool_kwargs['extra'] = {'biotoolsID': biotools_id}
    tool_kwargs.update(pop_websites_and_repo(meta_data['homepage'], meta_data['link'], meta_data['documentation']))
    tool_kwargs.update(_handle_credit(meta_data['credit']))
    return tool_kwargs

