# Copyright Sébastien Janel (2024- ) sebastien.janel@cnrs.fr
#
# This software is a computer program whose purpose is to analyze force curves
# recorded with an atomic force microscope.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use, modify
# and/ or redistribute the software under the terms of the CeCILL license as
# circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

"""Banner of the loading window: name, version, authors and logos of pyAF,
and a showcase of its results (images made with pyAF, in the colors of the
light and dark themes, crossfading)."""

import glob
import os

from qtpy import QtCore, QtGui, QtWidgets

from .. import consts
from ..gui_styles import theme_handler
from ..tools import misc_tools

# Height of the banner [logical pixels]: text on the left, showcase on the
# right in the rest of the width (the images are fitted in it, not cropped,
# with a margin)
BANNER_HEIGHT = 200
SHOWCASE_MARGIN = 10
RADIUS = 10
# Time an image is shown, and duration of the crossfade [ms]
SHOW_TIME = 4000
FADE_TIME = 900

# Colors of the banner (the background of the plots of the showcase)
COLORS = {"Light": {"background": "#ffffff", "border": "#d4d7dc",
                    "text": "#31363b", "secondary": "#6b7078"},
          "Dark": {"background": "#383C42", "border": "#4a4f56",
                   "text": "#eff0f1", "secondary": "#a6abb2"}}

# On two lines (not cut in the middle of a name)
AUTHORS = ("M. Popoff, S. Janel, V. Dupres, Y. Ciczora, S. Bovio,\n"
           "A. Dujardin, J. Lopez-Alonso, N. Duarte, F. Lafont")
# Color of the links, readable in both themes
LINK_COLORS = {"Light": "#1f6fd6", "Dark": "#7fb6f5"}


def showcase_images(theme):
    """Paths of the images of the showcase for a theme ("Light" or
    "Dark"), in their order."""
    folder = os.path.join(misc_tools.get_app_path(), "images", "showcase",
                          theme.lower())
    return sorted(glob.glob(os.path.join(folder, "*.png")))


class Showcase(QtWidgets.QWidget):
    """Images crossfading one into the next (no blank between them)."""

    def __init__(self, paths, background):
        super().__init__()
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Expanding)
        self.background = QtGui.QColor(background)
        self.pixmaps = []
        for path in paths:
            pixmap = QtGui.QPixmap(path)
            if not pixmap.isNull():
                self.pixmaps.append(pixmap)
        self.current = 0
        self._blend = 0.0

        self.animation = QtCore.QVariantAnimation(self)
        self.animation.setStartValue(0.0)
        self.animation.setEndValue(1.0)
        self.animation.setDuration(FADE_TIME)
        self.animation.setEasingCurve(QtCore.QEasingCurve.InOutQuad)
        self.animation.valueChanged.connect(self.set_blend)
        self.animation.finished.connect(self.next_image)

        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(SHOW_TIME)
        self.timer.timeout.connect(self.animation.start)
        if len(self.pixmaps) > 1:
            self.timer.start()

    def fitted_width(self, height):
        """Width needed to show the widest image at this height (with the
        margins)."""
        inner = height - 2 * SHOWCASE_MARGIN
        widths = [pixmap.width() * inner / pixmap.height()
                  for pixmap in self.pixmaps]
        return int(max(widths, default=inner)) + 2 * SHOWCASE_MARGIN

    def set_blend(self, value):
        self._blend = float(value)
        self.update()

    def next_image(self):
        self.current = (self.current + 1) % len(self.pixmaps)
        self._blend = 0.0
        self.update()

    def stop(self):
        """Stop the animation (window closed)."""
        self.timer.stop()
        self.animation.stop()

    def paintEvent(self, event):
        _ = event
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setRenderHint(QtGui.QPainter.SmoothPixmapTransform)
        # Rounded on the right, as the banner
        path = QtGui.QPainterPath()
        rect = QtCore.QRectF(self.rect())
        path.addRoundedRect(rect.adjusted(-RADIUS, 1, -1, -1), RADIUS - 1,
                            RADIUS - 1)
        painter.setClipPath(path)
        painter.fillRect(self.rect(), self.background)
        if not self.pixmaps:
            return
        self.draw_fitted(painter, self.pixmaps[self.current])
        if self._blend > 0:
            following = self.pixmaps[(self.current + 1) % len(self.pixmaps)]
            painter.setOpacity(self._blend)
            self.draw_fitted(painter, following)

    def draw_fitted(self, painter, pixmap):
        """Draw the whole image, as large as possible in the widget (with a
        margin), centered."""
        area = QtCore.QRectF(self.rect()).adjusted(
            SHOWCASE_MARGIN, SHOWCASE_MARGIN, -SHOWCASE_MARGIN,
            -SHOWCASE_MARGIN)
        scale = min(area.width() / pixmap.width(),
                    area.height() / pixmap.height())
        size = QtCore.QSizeF(pixmap.width() * scale, pixmap.height() * scale)
        target = QtCore.QRectF(QtCore.QPointF(0, 0), size)
        target.moveCenter(area.center())
        painter.drawPixmap(target, pixmap, QtCore.QRectF(pixmap.rect()))


class SplashBanner(QtWidgets.QFrame):
    """Banner of the loading window, in the colors of the current theme."""

    def __init__(self, width):
        super().__init__()
        theme = theme_handler.current_theme
        colors = COLORS.get(theme, COLORS["Dark"])
        self.setObjectName("pyaf_splash")
        self.setFixedSize(width, BANNER_HEIGHT)
        self.setStyleSheet(
            "QFrame#pyaf_splash { background-color: %s; border: 1px solid %s;"
            " border-radius: %dpx; }"
            "QFrame#pyaf_splash QLabel { background: transparent; "
            "border: none; color: %s; }"
            % (colors["background"], colors["border"], RADIUS,
               colors["text"]))

        # Name and version, with the icon
        icon = QtWidgets.QLabel()
        size = 44
        icon.setPixmap(QtGui.QIcon(
            os.path.join(misc_tools.get_app_path(), "images", "icon.svg"))
            .pixmap(size, size))
        title = QtWidgets.QLabel("pyAF")
        font = QtGui.QFont(self.font())
        font.setPointSize(26)
        font.setBold(True)
        title.setFont(font)
        version = QtWidgets.QLabel("Version " + consts.MAINVERSION)
        version.setStyleSheet("color: %s;" % colors["secondary"])
        names = QtWidgets.QVBoxLayout()
        names.setSpacing(0)
        names.addWidget(title)
        names.addWidget(version)
        # Logos, white text on the dark theme (images at 2x)
        logos = QtWidgets.QLabel()
        pixmap = QtGui.QPixmap(misc_tools.get_resource_path(
            "logos_" + ("dark" if theme == "Dark" else "light") + ".png"))
        pixmap.setDevicePixelRatio(2)
        logos.setPixmap(pixmap)
        logos.setFixedSize(pixmap.width() // 2, pixmap.height() // 2)

        head = QtWidgets.QHBoxLayout()
        head.setSpacing(12)
        head.addWidget(icon)
        head.addLayout(names)
        head.addStretch(1)

        tagline = QtWidgets.QLabel("An AFM Force Curve Analysis Software\n"
                                   "Dedicated to Biologists")
        font = QtGui.QFont(self.font())
        font.setPointSize(font.pointSize() + 1)
        tagline.setFont(font)
        authors = QtWidgets.QLabel(AUTHORS)
        font = QtGui.QFont(self.font())
        font.setItalic(True)
        font.setPointSize(max(font.pointSize() - 1, 8))
        authors.setFont(font)
        authors.setStyleSheet("color: %s;" % colors["secondary"])

        text = QtWidgets.QVBoxLayout()
        text.setContentsMargins(18, 12, 16, 12)
        text.setSpacing(4)
        text.addLayout(head)
        text.addWidget(tagline)
        text.addWidget(authors)
        text.addStretch(1)
        text.addWidget(logos, 0, QtCore.Qt.AlignLeft)

        self.showcase = Showcase(showcase_images(theme),
                                 colors["background"])
        # As wide as the widest image at the height of the banner: the text
        # and the showcase centered together in the banner
        self.showcase.setFixedWidth(self.showcase.fitted_width(BANNER_HEIGHT))

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addStretch(1)
        layout.addLayout(text)
        layout.addWidget(self.showcase)
        layout.addStretch(1)

    def stop(self):
        """Stop the showcase (window closed)."""
        self.showcase.stop()
