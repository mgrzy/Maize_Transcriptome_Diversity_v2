#!/usr/bin/env python3
"""Generate the presentation QR code; requires qrcode and pypng."""
from pathlib import Path

import qrcode
from qrcode.image.pure import PyPNGImage

PAPER_URL = "https://www.biorxiv.org/content/10.64898/2026.02.19.706772v1"
output = Path(__file__).resolve().parents[3] / "notes/images/biorxiv_paper_qr.png"
qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_M,
                   box_size=10, border=4)
qr.add_data(PAPER_URL)
qr.make(fit=True)
qr.make_image(image_factory=PyPNGImage).save(str(output))
print(f"Saved {output}: {PAPER_URL}")
