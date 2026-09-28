from pathlib import Path
import yaml

class Config:
    def __init__(self, config_file):
        config_file = Path(config_file).resolve()

        with open(config_file, "r", encoding="utf-8") as f:
            self.data = yaml.safe_load(f)

        self.config_dir = config_file.parent

    @property
    def periodic_table(self):
        path = Path(self.data["paths"]["periodic_table"])

        if not path.is_absolute():
            path = self.config_dir / path
        return path.resolve()
    