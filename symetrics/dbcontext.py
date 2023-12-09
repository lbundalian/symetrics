import sqlite3
from sqlite3 import Error
import logging

class DbContext:
    def __init__(self, db_file):
        self._db = db_file
        self._conn = None

    def __enter__(self):
        try:
            self._conn = sqlite3.connect(self._db)
            return self
        except Error as e:
            logging.error(f"Unable to connect to database {self._db}: {e}")
            raise

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._conn:
            self._conn.close()