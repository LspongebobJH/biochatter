from ._python_interpreter import evaluate_python_code

def run_codes(code: str, state: dict[str, object]):
    """
    Run codes

    Parameters
    -----------
    code : str
        A single valid code snippet as a string.
    state: dict[str, object]
        A dictionary of variables to be used in the code snippet. E.g. {'sc': sc, 'adata': adata}
    """
    
    try:
        result = str(evaluate_python_code(code, state=state)[0])
    except Exception as e:
        return f"ERROR: {str(e)}", e
    return result, None