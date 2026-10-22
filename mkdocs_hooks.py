from __future__ import annotations


def on_config(config):
    # mknotebooks builds these paths with os.path.join, which yields backslashes on Windows
    config["extra_css"] = [p.replace("\\", "/") for p in config["extra_css"]]
    return config
