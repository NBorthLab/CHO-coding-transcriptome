import snakemake

def parse_snakefile(rule, path):
    workflow = snakemake.Workflow(snakefile = path, rerun_triggers = "mtime")
    workflow.include(path)

    rule_object = workflow.get_rule(rule)

    inputs = {name: file for name, file in rule_object._input.items()}
    outputs = {name: file for name, file in rule_object._output.items()}
    params = {name: file for name, file in rule_object._params.items()}

    return {"input": inputs, "output": outputs}

