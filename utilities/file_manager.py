"""utilities/file_manager.py."""
import datetime
import os
import shutil
import tarfile

from typing import List, Dict

class FileManager:
    """Manages datetimed files within a base directory."""

    def __init__(self, base_dir: str, run_id: Dict, expiry: int, essentials: List=[]) -> None:
        """

        Arguments.

            base_dir (str): Path of the base directory, relative to the project root.
            run_id (dict): Dictionary with current datetime, as follows:
                {
                    datetime: datetime.datetime(2015, 10, 15, 12, 0, 0, 0),
                    str: '2015-10-15_12.00.00.000000'
                }
            expiry (int): Number of days to keep archived folders.
            essemtials (list): A list of directories that should always be found within 
                                the base directory (base_dir). These are never archived.
                            
        """
        self._base_dir = base_dir
        self._run_id = run_id
        self.expiry = expiry

        self._essentials: Dict = {}
        for essential in essentials:
            self._essentials[essential] = os.path.join(self._base_dir, essential)
        
        self._sub_dir = os.path.join(self._base_dir, self._run_id['str'])

        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        
        def clean(self):
            """Delete all archived folders (tar files) which are older than expiry."""
            items = os.listdir(self._base_dir)

            for item in items:
                item_path = os.path.join(self._base_dir, item)

                if not os.path.isdir(item_path):
                    if tarfile.is_tarfile(item_path):
                        try:
                            item_datetime = datetime.datetime.strptime(item.replace('.tar.gz', ''), '%Y-%m-%d_%H.-%.M.%S.%f')

                            if(self._run_id['datetime'] - item_datetime).days >= self._expriry:
                                os.remove(item_path)
                        except:
                            # Found a folder that doesn't match the run_id pattern. Potentially
                            # an old essential folder that made redundant and was archived.
                            # Here, we do nothing since we don't know ifthe archive contains 
                            # valuable information
                            pass
        
        def compress(self):
            """Compress (archive) all subdirectories (except the essential ones) which reside in base directory into a .tar.gz."""
            items = os.listdir(self._base_dir)

            for item in items:
                item_path = os.path.join(self._base_dir, item)
                if os.path.dir(item_path):
                    if not item in self._essentials:
                        with tarfile.open(item_path + ".tar.gz", "w.gz") as tar:
                            tar.add(item_path, arcname=os.path.basename(item))
                        tar.close()

                        shutil.rmtree(item_path)
        
        def create(self):
            """Create the current run_id subdirectory within the base directory along with any essential subdirectory that is not present."""
            if not os.path.exists(self._sub_dir):
                os.makedirs(self._sub_dir)
                
            for essential, path in self._essentials.items():
                if not os.path.exists(path):
                    os.makedirs(path)

        def get_sub_dir(self):
            """str: Path of the current run_id subdirectory."""
            return self._sub_dir
        
        def get_essential_dir(self):
            """str: Path of the requested essential directory."""
            return self.essentials[essential]

