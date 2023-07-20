
import hashlib
import logging
import os
import subprocess as sp

log_format = "%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s"

def init_log(log_file_path=None, rootloglevel=logging.DEBUG):
    default_formatter = logging.Formatter(log_format)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(default_formatter)
    root = logging.getLogger()

    if log_file_path:
        fh = logging.FileHandler(log_file_path)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(default_formatter)
        root.addHandler(fh)

    root.addHandler(console_handler)
    root.setLevel(rootloglevel)

_log = logging.getLogger(__name__)

def mkdir(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    
    if not os.path.exists(dirpath):
        err = f'"{dirpath}" could not be created'
        _log.error(err)
        raise FileNotFoundError(err)


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.__dict__)

def execute(cmd, retcodes=[0], stdout=sys.stdout, stderr=sys.stderr, exec_mode=None):
    cmd = cmd.strip()
    exec_mode = exec_mode if exec_mode else os.environ.get("TP_EXEC_MODE", "raw")
    if (exec_mode == "print"):
        stdout.write(cmd)
        return 0
    try:
        _log.debug(cmd)
        process = sp.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
        process.communicate()
        return process.returncode

    except CalledProcessError as ex:
        if ex.returncode not in retcodes:
            _log.critical(ex, exc_info=True)
            raise
    return ex.returncode

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def md5_equal(fname, known_md5):
    match = md5(fname) == known_md5
    if not match:
        _log.error("error in %s checksum" % fname)
    return match


def download_file(complete_url, target, ovewrite=False, retries=3,timeout=20):
    if not target.strip():
        target = "./"
    if not os.path.exists(os.path.dirname(os.path.abspath(target))):
        raise FileNotFoundError("%s does not exists" % os.path.dirname(target))
    if os.path.exists(target) and not ovewrite:
        raise OvewriteFileException("%s already exists" % target)

    execute(f'wget  --timeout={timeout} --tries={retries} -O "{target}" "{complete_url}"')
