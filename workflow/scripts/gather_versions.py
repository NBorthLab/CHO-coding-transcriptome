import glob
import yaml
from yaml import CLoader, CDumper
import sys

# %%

def get_version_yamls() -> list[str]:
    return glob.glob("logs/**/*_versions.yaml", recursive = True)


def cat_yamls(yamls: list[str]) -> dict[str, str]:

    all_versions = {}

    for file in yamls:
        with open(file) as f:
            file_versions: dict(str, str) = yaml.load(f, Loader = CLoader)
        all_versions.update(file_versions)

    return all_versions

def print_yaml(versions: dict[str, str]) -> None:
    with open("logs/versions.yaml", "w") as f:
        yaml.dump(versions, f, CDumper)

# %%

def main() -> None:
    yaml_files: list[str] = get_version_yamls()
    all_versions = cat_yamls(yaml_files)
    print_yaml(all_versions)

# %%

if __name__ == "__main__":
    main()

