
def get_dict_from_list(dict_list, target_key_name, value):
    """
    Find a dictionary in a list dictionaries that has a particular value for a key.
    :param dict_list:
    :param target_key_name:
    :param value:
    :param parent_key: if the key you are looking for is nested in a parent field.
    :return:
    """
    items_found = 0
    target_index = -1
    for index, dict_item in enumerate(dict_list):
        # skip dict if key doesn't exist.
        curr_dict_item = dict_item
        item_value = curr_dict_item.get(target_key_name)
        if item_value:  # key exists in dict and has a value.
            if curr_dict_item[target_key_name] == value:
                my_dict = dict_item
                target_index = index
                items_found += 1
    if items_found == 0:
        raise ValueError(f"No dictionary was found for key: {target_key_name}, value: {value}")
    elif items_found > 1:
        raise ValueError(f"Found {items_found} dictionaries with key: {target_key_name}, value: {value}. Value is not unique")
    elif items_found == 1:
        return my_dict, target_index
    else:
        raise NotImplementedError("Should never hit this.")