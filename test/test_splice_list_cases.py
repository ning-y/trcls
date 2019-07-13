import os, sys

sys.path.append(os.path.join(sys.path[0], '../src'))

from splice_list import Junction

exon_list = [
        (1, 2), (412, 532), (9999999999999998, 9999999999999999)
]

exons_list = [
        (
            'TEST_SPLICE_VARIANT',
            ((1, 2), (412, 532), (9999999999999998, 9999999999999999))
        )
]

overlapping_exons_list = [
        (
            'TEST_SPLICE_VARIANT',
            ((1, 2), (2, 3))
        ),
        (
            'TEST_SPLICE_VARIANT',
            ((2, 3), (1, 2))
        ),
        (
            'TEST_SPLICE_VARIANT',
            ((200, 300), (100, 250))
        )
]

junction_list = [
        (1, Junction.TYPE_START), (2, Junction.TYPE_END),
        (412, Junction.TYPE_START), (532, Junction.TYPE_END),
        (9999999999999998, Junction.TYPE_START),
        (9999999999999999, Junction.TYPE_END)
]

matching_exons = [
        [
            [(1,10), (11, 20), (31, 40)],
            [(1,10), (11, 20), (31, 40)],
            0
        ],
        [
            [(1,10), (11, 20), (31, 40)],
            [(2,9), (10, 21), (30, 39)],
            1
        ],
        [
            [(1,10), (11, 20), (31, 40)],
            [(6,15), (16, 20), (33, 44)],
            5
        ],
        [
            [(1,10), (11, 20), (31, 40)],
            [(11, 20)],
            0
        ]
]

non_matching_exons = [
        [
            [(1,10), (31, 40)],
            [(1,10), (11, 20), (31, 40)],
            0
        ],
        [
            [(1,10), (11, 20), (31, 40)],
            [(1,10), (11, 21), (31, 40)],
            0
        ],
        [
            [(1,10), (11, 20), (31, 40)],
            [(1,5), (11, 21), (31, 40)],
            0
        ],
        [
            [(1,10), (13, 20), (31, 40)],
            [(2,12), (13, 21), (31, 40)],
            1
        ],
]
