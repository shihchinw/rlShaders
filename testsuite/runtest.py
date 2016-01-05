#!/usr/bin/env
#
# rlShader test utility.
# Please execute following command to see the usages.
# $> python runtest.py -h
import argparse
import csv
import glob
import os
import platform
import re
import subprocess
import sys


def enum(**enums):
    return type('Enum', (), enums)

_TASK = enum(MKDIR=1, RENDER=2, LIST=3, DISPLAY=4)
_COLOR_CODE = enum(STATUS='\033[92m', WARNING='\033[93m', ERROR='\033[91m', END='\033[0m')


def parse_args(args):
    epilog = r"""Examples:

    mkdir --sn 1..3
    => Create folder structure for test cases 1, 2, 3.

    render --ai C:\solidangle\mtoadeploy\2015 -l extra-shader-path --as-ref --sn 1 2
    => Render reference images for test cases 1 & 2.

    render --ai C:\solidangle\mtoadeploy\2015 -l extra-shader-path --sn 1..3 7
    => Render and diff images for test cases 1, 2, 3 & 7.

    render --ai C:\solidangle\mtoadeploy\2015 -l extra-shader-path
    => Render and diff images for ALL test cases.

    list
    => List test case description.

    display -l extra-shader-path --sn 1 3 5
    => Use iv from OIIO to view images of test cases 1, 3, 5.
    """

    parser = argparse.ArgumentParser(description='rlShader test utility.', epilog=epilog,
                                     formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers()
    parser_mkdir = subparsers.add_parser('mkdir', help='Create test folders')
    parser_mkdir.set_defaults(task_type=_TASK.MKDIR)

    parser_render = subparsers.add_parser('render', help='Render test cases')
    parser_render.add_argument('--ai', type=str, dest="arnold_path", required=True,
                               help='Arnold path')
    parser_render.add_argument('-l', type=str, dest="shader_path", help='Shader search path')
    parser_render.add_argument('--as-ref', dest="as_ref", action='store_true',
                               help='Render image as reference')
    parser_render.add_argument('-v', dest="verbose", action='store_true', help='Verbose')
    parser_render.set_defaults(task_type=_TASK.RENDER)

    parser_list = subparsers.add_parser('list', help='List test cases')
    parser_list.set_defaults(task_type=_TASK.LIST)
    parser_display = subparsers.add_parser('display', help='Display reference/test images')
    parser_display.set_defaults(task_type=_TASK.DISPLAY)

    required = {'mkdir': True, 'render': False, 'list': False, 'display': True}
    for k, p in subparsers.choices.iteritems():
        p.add_argument("--sn", type=str, dest="serial", nargs='+',
                       required=required[k], help='Test case number(s)')

    if args:
        return parser.parse_args(args)

    print parser.print_help()


def _prompt(msg, color_code=None):
    if not color_code:
        color_code = _COLOR_CODE.STATUS
    print '%s%s%s' % (color_code, msg, _COLOR_CODE.END)


def create_folder(case_id, root):
    """Create folder structure for an unit test.

    [root]
      |--[data]
      |--[ref]
      |--README
    """
    testname = str(case_id).zfill(4)
    testpath = os.path.join(root, testname)

    if os.path.exists(testpath):
        _prompt('"%s" already exists' % testpath, _COLOR_CODE.WARNING)
        return

    os.mkdir(testpath)
    os.mkdir(os.path.join(testpath, 'data'))
    os.mkdir(os.path.join(testpath, 'ref'))

    with open(os.path.join(testpath, 'README'), 'w') as f:
        f.write('TODO: test description')
    _prompt('Create "%s"' % testpath)


def expand_serial_no(serial_strings):
    """Retrieve serial numbers from a list of strings.

    Ex. ["2..5", "8"] => [2, 3, 4, 5, 8]
    """
    case_ids = []
    for serial in serial_strings:
        result = re.match('(\d+)[.]{2}(\d+)', serial)
        if result:
            start, end = map(int, result.groups())
            case_ids.extend(range(start, end + 1))
        else:
            case_ids.append(int(serial))

    return case_ids


class Test(object):

    @staticmethod
    def run_tests(testdir, case_ids, as_ref=False, env=None, verbose=False):
        report_file_path = os.path.join(testdir, 'report.csv')
        with open(report_file_path, 'wb') as csvfile:
            report_writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)

            for case_id in case_ids:
                t = Test(testdir, case_id, env, verbose)
                result = t.render(as_ref)
                if not as_ref:
                    result = t.compare()

                status = ['FAIL', 'OK'][result]
                color_code = [_COLOR_CODE.ERROR, _COLOR_CODE.STATUS][result]
                _prompt('%s (%s)' % (t, status), color_code)
                report_writer.writerow([case_id, t.get_description(), status])

    @staticmethod
    def list_tests(testdir, case_ids):
        print '-' * 80
        print 'Tests in "%s"' % testdir
        print '-' * 80

        for case_id in case_ids:
            print Test(testdir, case_id)

    @staticmethod
    def display_tests(testdir, case_ids, env=None):
        for case_id in case_ids:
            Test(testdir, case_id, env).display()

    def __init__(self, root, case_id, env=None, verbose=False):
        self.case_id = str(case_id).zfill(4)
        self.root = os.path.join(root, self.case_id)
        self.data_dir = os.path.join(self.root, 'data')
        self.ref_dir = os.path.join(self.root, 'ref')
        self.env = env
        self.testfile = self._get_file()
        self.verbose = verbose

    def __str__(self):
        return '[%s] %s' % (self.case_id, self.get_description())

    def _get_file(self):
        testfile = glob.glob('%s/data/*.ass' % self.root)
        if testfile and len(testfile) == 1:
            return os.path.normpath(testfile[0])

    def _prompt(self, msg, color_code):
        _prompt('[%s] %s' % (self.case_id, msg), color_code)

    def _prompt_status(self, msg):
        if self.verbose:
            self._prompt(msg, _COLOR_CODE.STATUS)

    def _prompt_warning(self, msg):
        if self.verbose:
            self._prompt(msg, _COLOR_CODE.WARNING)

    def _prompt_error(self, msg):
        self._prompt(msg, _COLOR_CODE.ERROR)

    def get_description(self):
        descfile = os.path.join(self.root, 'README')
        with open(descfile, 'r') as f:
            return f.readline().strip()

    def render(self, as_ref=False):
        if not self.testfile:
            self._prompt_error('Test file not found!')
            return False

        texture_path = self.data_dir + ';' + self.env['TEXPATH']
        render_args = ['kick -dw -dp -set options.texture_searchpath "%s"' % texture_path]

        if 'SHADERPATH' in self.env:
            render_args.append('-l "%s"' % self.env['SHADERPATH'])

        token = ['test', 'ref'][as_ref]
        imgpath = os.path.join(self.ref_dir, '%s.exr' % token)
        render_args.append(' -i "%s" -o "%s"' % (self.testfile, imgpath))
        cmd = ' '.join(render_args)

        self._prompt_status('Rendering...')
        logfile = os.path.join(self.ref_dir, '%s.log' % token)
        with open(logfile, 'w') as f:
            retcode = subprocess.call(cmd, shell=True, env=self.env, stdout=f)
            if retcode != 0:
                self._prompt_error('Failed to render with flags: "%s"' % cmd)
                print >>sys.stderr
                return False

        return True

    def compare(self, max_rms_error=0.005):
        """Compare render result with idiff from OpenImageIO.

        Returns:
            bool: False if the test image is not found or the RMS error is less than max_rms_error.
        """

        try:
            ref_img, test_img = glob.glob('%s/ref/*.exr' % self.root)
        except:
            self._prompt_warning('Missing file for comparison...')
            return False

        self._prompt_status('Comparing...')

        cmd = 'idiff "%s" "%s"' % (os.path.normpath(ref_img), os.path.normpath(test_img))
        try:
            subprocess.check_output(cmd, shell=True, env=self.env)
            return True
        except subprocess.CalledProcessError as inst:
            if self.verbose:
                print inst.output
            m = re.search('\sRMS error = (\d+[.]\d+)', inst.output)
            rms_error = float(m.group(1))
            return rms_error < max_rms_error

    def display(self):
        images = glob.glob('%s/ref/*.exr' % self.root)

        if not images:
            self._prompt_warning('Missing image files')
            return

        cmd = 'iv %s' % (' '.join(images))
        subprocess.call(cmd, shell=True, env=self.env)


def main(args, test_root, suite='mtoa'):
    os.chdir(test_root)

    suite_dir = os.path.join(test_root, suite)

    if not args.serial:
        case_ids = filter(lambda x: re.match('\d{4}', x), os.listdir(suite_dir))
    else:
        case_ids = expand_serial_no(args.serial)

    if args.task_type is _TASK.MKDIR:
        map(lambda x: create_folder(x, suite_dir), case_ids)
    elif args.task_type is _TASK.LIST:
        Test.list_tests(suite_dir, case_ids)
    elif args.task_type is _TASK.DISPLAY:
        test_environ = os.environ.copy()

        if (platform.system() == 'Windows'):
            test_environ['PATH'] = '%s;%s' % (os.path.join(test_root, 'bin'), test_environ['PATH'])
        Test.display_tests(suite_dir, case_ids, test_environ)
    else:
        shader_path = os.path.join(args.arnold_path, 'shaders')
        if args.shader_path:
            shader_path += ';' + args.shader_path.rstrip('/\\')

        texture_path = os.path.join(test_root, "data")
        test_environ = {'SHADERPATH': shader_path, 'TEXPATH': texture_path}
        test_environ.update(os.environ)

        if (platform.system() == 'Windows'):
            binary_dirs = [os.path.join(args.arnold_path, 'bin'),
                           os.path.join(test_root, 'bin')]

            test_environ['PATH'] = '%s;%s' % (';'.join(binary_dirs), test_environ['PATH'])
        Test.run_tests(suite_dir, case_ids, args.as_ref, test_environ, args.verbose)

if __name__ == '__main__':
    test_root = os.path.dirname(os.path.realpath(__file__))
    main(parse_args(sys.argv[1:]), test_root)
