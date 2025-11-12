"""
Components of the Graphical user interface
    (Personalized PyQt objects)

adapted from David Patsch's AlphaDock

components.py
"""

from . import styles, config

import webbrowser

from PyQt5.QtWidgets import (QHBoxLayout, QPushButton, QComboBox, QWidget, 
                             QLayout, QMenu, QAction, QCheckBox, QLineEdit, 
                             QWidgetAction)
from PyQt5.QtCore import pyqtSignal, Qt, QPoint, QEvent

from functools import partial
import re

class ComboBox(QComboBox):
    """Custom QComboBox emitting a signal by hovering over"""
    entered = pyqtSignal()

    def enterEvent(self, event: QEvent) -> None:
        super().enterEvent(event)
        self.entered.emit()


class RightPushButton(QPushButton):
    """
    QPushbutton with additional function on right click
    adapted from stackoverflow.com/questions/44264157
    """
    def __init__(self, parent: QLayout) -> None:
        super(RightPushButton, self).__init__()
        self.parent = parent

    def mousePressEvent(self, QMouseEvent: QEvent) -> None:
        """If rightclicked, call rightClickedEvent function from InteractionGUI"""
        if QMouseEvent.button() == Qt.RightButton:
            self.parent.parent.rightClickedEvent()


class QAMenuHosts(QMenu):
    """
    QMenu adapted to add and handle QAction widgets for the different hosts
    """
    def __init__(self, hosts: dict, parent: QWidget) -> None:
        super().__init__("hosts", parent)

        # initialize general variables
        self.hosts = hosts
        self.parent = parent

        # create QAction widgets for the different hosts
        for t, (k, v) in enumerate(self.hosts.items()):
            action = QAction(v["alias"], self, triggered=partial(
                self.change_host_and_update_ticks, k, v, t))
            action.setCheckable(True)
            self.addAction(action)
            if t == 0:
                self.change_host_and_update_ticks(k, v, t)

    def change_host_and_update_ticks(self, new_host: str, configs: dict, idx: int) -> None:
        """
        change host and update ticks
        
        :param str new_host:                        hostname
        :param dict[str[Union[str, int]]] configs:  host parameters (alias, num_cpu, port)
        :param int idx:                             index of QAction widget (linked to host)
        """
        self.parent.parent.change_host(new_host, configs)
        for act in self.actions():
            act.setChecked(False)
        self.actions()[idx].setChecked(True)


class QAMenuIntraChains(QMenu):
    """
    QMenu adapted to add and handle QAction widgets for intra chain interactions
    """
    def __init__(self, parent: QWidget) -> None:
        super().__init__("Intra chain interactions", parent)
        self.setToolTip("If enabled: intra-chain interactions are analyzed\n"
                        "By default only pdb files are analyzed")

        # initialize general variables
        self.parent = parent

        # create chain input
        self.chainCbox = QCheckBox()
        self.chainCbox.clicked.connect(self.setChain)

        self.chainLEdit = QLineEdit()
        self.chainLEdit.setMinimumWidth(50)
        self.chainLEdit.textChanged.connect(self.change_chain)
        self.chainLEdit.setToolTip("Specify chain e.g., A")
        self.chainLEdit.setEnabled(self.chainCbox.isChecked())
        
        chainLayout = QHBoxLayout()
        chainLayout.addWidget(self.chainCbox)
        chainLayout.addWidget(self.chainLEdit)
        
        chainWidget = QWidget()
        chainWidget.setLayout(chainLayout)
        
        chainWAction = QWidgetAction(self)
        chainWAction.setDefaultWidget(chainWidget)
        self.addAction(chainWAction)

    def setChain(self) -> None:
        """
        set chain
        """
        self.chainLEdit.setEnabled(self.chainCbox.isChecked())
        if not self.chainCbox.isChecked():
            self.parent.parent.API.intraChain = "0"
        else:
            if not self.parent.parent.pymolBtn.isChecked():
                self.parent.parent.pdbBtn.setChecked(True)

    def change_chain(self) -> None:
        """
        update chain
        """
        chain = self.chainLEdit.text().upper()

        # input validation & exception handling
        if len(chain) != 1 or not re.match("[A-Z]", chain):
            print(f"{self.chainLEdit.text()} is not a valid chain")
            chain = "0"
        
        self.parent.parent.API.intraChain = chain
        
        
class MyBar(QWidget):
    """QWidget to replace the removed title bar and buttons"""

    def __init__(self, parent: QLayout) -> None:
        super(MyBar, self).__init__()
        
        # initialize general variables
        self.setStyleSheet(styles.David)      
        self.parent = parent
        
        # general settings
        layout = QHBoxLayout()
        layout.setContentsMargins(0,0,0,0)
        btn_size = 30

        # left hand side
        leftWidget = QWidget()
        leftLayout = QHBoxLayout()
        leftLayout.setContentsMargins(0,0,0,0)

        # logo
        ccbiologo = RightPushButton(self)
        ccbiologo.setStyleSheet(styles.ccbiologo)
        ccbiologo.setFixedSize(btn_size,btn_size)

        # settings
        btn_settings = QPushButton()
        btn_settings.setToolTip('Settings (Ctrl+I)')
        btn_settings.setShortcut("Ctrl+I")
        btn_settings.setFixedSize(btn_size,btn_size)
        btn_settings.setStyleSheet(styles.SettingsButton)

        self.actionshowDist = QAction("interaction lengths", self)
        self.actionshowDist.setShortcut("Ctrl+D")
        self.actionshowDist.setCheckable(True)

        self.actionclean_up = QAction("reset everything", self)
        self.actionclean_up.setShortcut("Ctrl+R")

        submenu_hosts = QAMenuHosts(config.HOSTS, self)
        self.submenu_chains = QAMenuIntraChains(self)

        menu_setting = QMenu()
        menu_setting.setToolTipsVisible(True)
        menu_setting.setStyleSheet(styles.David)
        menu_setting.addAction(self.actionclean_up)
        menu_setting.addAction(self.actionshowDist)
        menu_setting.addMenu(submenu_hosts)
        menu_setting.addMenu(self.submenu_chains)
        btn_settings.setMenu(menu_setting)

        # help
        btn_help = QPushButton()
        btn_help.clicked.connect(lambda: webbrowser.open(config.HELP_PATH))
        btn_help.setToolTip('Help (Ctrl+H)')
        btn_help.setFixedSize(btn_size, btn_size)
        btn_help.setStyleSheet(styles.MapButton)
        btn_help.setShortcut("Ctrl+H")

        leftLayout.addWidget(ccbiologo)
        leftLayout.addWidget(btn_settings)
        leftLayout.addWidget(btn_help)
        leftWidget.setLayout(leftLayout)
        layout.addWidget(leftWidget, alignment=Qt.AlignLeft)

        # right hand side
        rightWidget = QWidget()
        rightLayout = QHBoxLayout()
        rightLayout.setContentsMargins(0,0,0,0)

        # minimize window
        btn_min = QPushButton()
        btn_min.clicked.connect(lambda: self.parent.showMinimized())
        btn_min.setFixedSize(btn_size, btn_size)
        btn_min.setStyleSheet(styles.MinButton)

        # close window
        btn_close = QPushButton()
        btn_close.clicked.connect(lambda: self.parent.close())
        btn_close.setFixedSize(btn_size,btn_size)
        btn_close.setStyleSheet(styles.CloseButton)

        rightLayout.addWidget(btn_min)
        rightLayout.addWidget(btn_close)
        rightWidget.setLayout(rightLayout)
        layout.addWidget(rightWidget, alignment=Qt.AlignRight)

        self.setLayout(layout)

        # ability to move the window
        # adapted from from stackoverflow.com/questions/44241612
        self.start = QPoint(0, 0)
        self.pressing = False

    def mousePressEvent(self, event: QEvent) -> None:
        """
        ability to move the window
        adapted from from stackoverflow.com/questions/44241612
        """
        self.start = self.mapToGlobal(event.pos())
        self.pressing = True

    def mouseMoveEvent(self, event: QEvent) -> None:
        """
        ability to move the window
        adapted from from stackoverflow.com/questions/44241612
        """
        if self.pressing:
            self.end = self.mapToGlobal(event.pos())
            self.movement = self.end-self.start
            self.parent.setGeometry(self.mapToGlobal(self.movement).x()-10,
                                self.mapToGlobal(self.movement).y()-10,
                                self.parent.width(),
                                self.parent.height())
            self.start = self.end
