from .setup import TestCase


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_me(self):
        pass
