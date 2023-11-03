import sys

sys.path.append('release')

import evogen
import unittest

class TestEVO( unittest.TestCase ):

    def testSetup( self ):

        res = evogen.add(1, 2)

        self.assertEqual( res, 3 )

if __name__ == '__main__':
    unittest.main()
    