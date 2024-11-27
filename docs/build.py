from pathlib import Path

here = Path(__file__).parent
import build_scripts.build_art_api
import build_scripts.build_artcore_api
# Touch the addons-api.md file to trigger a rebuild of the addons API page
(here / "src" / "content" / "modified-api.md").touch()
