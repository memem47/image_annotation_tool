from main import _scale_to_size
def test_scale():
    assert _scale_to_size((100, 100), 0.25) == (25, 25)
    assert _scale_to_size((200, 50), 4.0) == (800, 800)