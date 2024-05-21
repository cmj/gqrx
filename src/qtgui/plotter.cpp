/* -*- c++ -*- */
/* + + +   This Software is released under the "Simplified BSD License"  + + +
 * Copyright 2010 Moe Wheatley. All rights reserved.
 * Copyright 2011-2013 Alexandru Csete OZ9AEC
 *
 * Redistribution and use in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice, this list of
 *       conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list
 *       of conditions and the following disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY Moe Wheatley ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Moe Wheatley OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those of the
 * authors and should not be interpreted as representing official policies, either expressed
 * or implied, of Moe Wheatley.
 */
#include <cmath>
#include <QColor>
#include <QDateTime>
#include <QDebug>
#include <QFont>
#include <QPainter>
#include <QtGlobal>
#include <QToolTip>
#include "plotter.h"
#include "bandplan.h"
#include "bookmarks.h"
#include "dxc_spots.h"
#include <volk/volk.h>

Q_LOGGING_CATEGORY(plotter, "plotter")

#define CUR_CUT_DELTA         5     // cursor capture delta in pixels
#define CLICK_FREQ_RESOLUTION 100   // frequency rounding for set via click
#define VDIV_DELTA            30

#define FFT_MIN_DB     -160.f
#define FFT_MAX_DB      30.f
#define FFT_MIN_DB_RANGE 2.f

#define FILTER_WIDTH_MIN_HZ 200

// Colors of type QRgb in 0xAARRGGBB format (unsigned int)
#define PLOTTER_BGD_COLOR           0xFF1F1D1D
#define PLOTTER_GRID_COLOR          0x80606060
#define PLOTTER_TEXT_COLOR          0xFFDADADA
#define PLOTTER_CENTER_LINE_COLOR   0x80CCDDFF
#define PLOTTER_FILTER_LINE_COLOR   0xB0FF6060
#define PLOTTER_FILTER_BOX_COLOR    0x28FFFFFF
#define PLOTTER_MARKER_COLOR        0XB080FF80
// FIXME: Should cache the QColors also

#define HOR_MARGIN 5
#define VER_MARGIN 5

static inline bool val_is_out_of_range(float val, float min, float max)
{
    return (val < min || val > max);
}

static inline bool out_of_range(float min, float max)
{
    return (val_is_out_of_range(min, FFT_MIN_DB, FFT_MAX_DB) ||
            val_is_out_of_range(max, FFT_MIN_DB, FFT_MAX_DB) ||
            max < min + FFT_MIN_DB_RANGE);
}

#define STATUS_TIP \
    "Click, drag or scroll on spectrum to tune. " \
    "Drag and scroll X and Y axes for pan and zoom. " \
    "Drag filter edges to adjust filter."

CPlotter::CPlotter(QWidget *parent) : QFrame(parent)
{
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setFocusPolicy(Qt::StrongFocus);
    setAttribute(Qt::WA_PaintOnScreen,false);
    setAutoFillBackground(false);
    setAttribute(Qt::WA_OpaquePaintEvent, false);
    setAttribute(Qt::WA_NoSystemBackground, true);
    setMouseTracking(true);

    m_TooltipsEnabled = false;
    setStatusTip(tr(STATUS_TIP));
    setWfColormap("gqrx");

    m_MaxHoldActive = false;
    m_MaxHoldValid = false;
    m_MinHoldActive = false;
    m_MinHoldValid = false;
    m_IIRValid = false;
    m_histIIRValid = false;
    m_alpha = 1.0;
    m_histMaxIIR = std::numeric_limits<float>::min();

    m_FftCenter = 0;
    m_CenterFreq = 144500000;
    m_DemodCenterFreq = 144500000;
    m_DemodHiCutFreq = 5000;
    m_DemodLowCutFreq = -5000;

    m_FLowCmin = -25000;
    m_FLowCmax = -1000;
    m_FHiCmin = 1000;
    m_FHiCmax = 25000;
    m_symetric = true;

    m_ClickResolution = CLICK_FREQ_RESOLUTION;
    m_FilterClickResolution = CLICK_FREQ_RESOLUTION;
    m_CursorCaptureDelta = CUR_CUT_DELTA;
    m_WaterfallMode = WATERFALL_MODE_MAX;
    m_PlotMode = PLOT_MODE_MAX;
    m_PlotScale = PLOT_SCALE_DBFS;
    m_PlotPerHz = false;

    m_FilterBoxEnabled = true;
    m_CenterLineEnabled = true;
    m_MarkersEnabled = false;
    m_BandPlanEnabled = true;
    m_BookmarksEnabled = true;
    m_InvertScrolling = false;
    m_DXCSpotsEnabled = true;

    m_Span = 96000;
    m_SampleFreq = 96000;

    m_HorDivs = 12;
    m_VerDivs = 6;
    m_PandMaxdB = m_WfMaxdB = 0.f;
    m_PandMindB = m_WfMindB = FFT_MAX_DB;

    m_CumWheelDelta = 0;
    m_FreqUnits = 1000000;
    m_CursorCaptured = NOCAP;
    m_Running = false;
    m_DrawOverlay = true;
    m_2DPixmap = QPixmap();
    m_OverlayPixmap = QPixmap();
    m_WaterfallImage = QImage();
    m_Size = QSize(0,0);
    m_GrabPosition = 0;
    m_Percent2DScreen = 35;	//percent of screen used for 2D display
    m_VdivDelta = VDIV_DELTA;
    m_BandPlanHeight = 0.0;

    m_FreqDigits = 6;

    m_Peaks = QMap<int,qreal>();
    enablePeakDetect(false);

    setFftPlotColor(QColor(0xFF,0xFF,0xFF,0xFF));
    enableFftFill(false);

    // always update waterfall
    tlast_wf_ms = 0;
    tlast_plot_drawn_ms = 0;
    tlast_wf_drawn_ms = 0;
    wf_valid_since_ms = 0;
    msec_per_wfline = 0;
    tlast_peaks_ms = 0;
    wf_epoch = 0;
    wf_count = 0;
    wf_avg_count = 0;
    wf_span = 0;
    fft_rate = 15;
}

CPlotter::~CPlotter()
= default;

QSize CPlotter::minimumSizeHint() const
{
    return {50, 50};
}

QSize CPlotter::sizeHint() const
{
    return {180, 180};
}

void CPlotter::mouseMoveEvent(QMouseEvent* event)
{
    QPoint pt = event->pos();

    int w = m_OverlayPixmap.width();
    int h = m_OverlayPixmap.height();
    int px = qRound((qreal)pt.x() * m_DPR);
    int py = qRound((qreal)pt.y() * m_DPR);
    QPoint ppos = QPoint(px, py);

    /* mouse enter / mouse leave events */
    if (py < h)
    {
        //is in Overlay bitmap region
        if (event->buttons() == Qt::NoButton)
        {
            bool onTag = false;
            if(py < 15 * 10) // FIXME
            {
                if(m_BookmarksEnabled || m_DXCSpotsEnabled)
                {
                    for(int i = 0; i < m_Taglist.size() && !onTag; i++)
                    {
                        if (m_Taglist[i].first.contains(ppos))
                            onTag = true;
                    }
                }
            }
            // if no mouse button monitor grab regions and change cursor icon
            if (onTag)
            {
                setCursor(QCursor(Qt::PointingHandCursor));
                m_CursorCaptured = TAG;
            }
            else if (isPointCloseTo(px, m_YAxisWidth/2, m_YAxisWidth/2))
            {
                if (YAXIS != m_CursorCaptured)
                    setCursor(QCursor(Qt::OpenHandCursor));
                m_CursorCaptured = YAXIS;
                if (m_TooltipsEnabled)
                    QToolTip::hideText();
            }
            else if (isPointCloseTo(py, m_XAxisYCenter, m_CursorCaptureDelta+20))
            {
                if (XAXIS != m_CursorCaptured)
                    setCursor(QCursor(Qt::OpenHandCursor));
                m_CursorCaptured = XAXIS;
                if (m_TooltipsEnabled)
                    QToolTip::hideText();
            }
            else if (isPointCloseTo(px, m_DemodFreqX, m_CursorCaptureDelta))
            {
                // in move demod box center frequency region
                if (CENTER != m_CursorCaptured)
                    setCursor(QCursor(Qt::SizeHorCursor));
                m_CursorCaptured = CENTER;
                if (m_TooltipsEnabled)
                    showToolTip(event, QString("Demod: %1 kHz").arg(m_DemodCenterFreq/1.e3, 0, 'f', 3));
            }
            else if (isPointCloseTo(px, m_DemodHiCutFreqX, m_CursorCaptureDelta))
            {
                // in move demod hicut region
                if (RIGHT != m_CursorCaptured)
                    setCursor(QCursor(Qt::SizeFDiagCursor));
                m_CursorCaptured = RIGHT;
                if (m_TooltipsEnabled)
                    showToolTip(event, QString("High cut: %1 Hz").arg(m_DemodHiCutFreq));
            }
            else if (isPointCloseTo(px, m_DemodLowCutFreqX, m_CursorCaptureDelta))
            {
                // in move demod lowcut region
                if (LEFT != m_CursorCaptured)
                    setCursor(QCursor(Qt::SizeBDiagCursor));
                m_CursorCaptured = LEFT;
                if (m_TooltipsEnabled)
                    showToolTip(event, QString("Low cut: %1 Hz").arg(m_DemodLowCutFreq));
            }
            else if (m_MarkersEnabled && isPointCloseTo(px, m_MarkerAX, m_CursorCaptureDelta))
            {
                if (MARKER_A != m_CursorCaptured && m_MarkerFreqA != MARKER_OFF)
                    setCursor(QCursor(Qt::OpenHandCursor));
                m_CursorCaptured = MARKER_A;
                if (m_TooltipsEnabled)
                    showToolTip(event, QString("Marker A: %1 kHz").arg(m_MarkerFreqA/1.e3, 0, 'f', 3));
            }
            else if (m_MarkersEnabled && isPointCloseTo(px, m_MarkerBX, m_CursorCaptureDelta))
            {
                if (MARKER_B != m_CursorCaptured && m_MarkerFreqB != MARKER_OFF)
                    setCursor(QCursor(Qt::OpenHandCursor));
                m_CursorCaptured = MARKER_B;
                if (m_TooltipsEnabled)
                    showToolTip(event, QString("Marker B: %1 kHz").arg(m_MarkerFreqB/1.e3, 0, 'f', 3));
            }
            else
            {	//if not near any grab boundaries
                if (NOCAP != m_CursorCaptured)
                {
                    setCursor(QCursor(Qt::ArrowCursor));
                    m_CursorCaptured = NOCAP;
                }
                if (m_TooltipsEnabled)
                {
                    QString toolTipText;
                    qint64 hoverFrequency = freqFromX(px);
                    toolTipText = QString("%1 kHz\nΔ %2 kHz")
                                          .arg(hoverFrequency/1.e3, 0, 'f', 3)
                                          .arg(locale().toString((hoverFrequency - m_DemodCenterFreq)/1.e3, 'f', 3));

                    QFontMetricsF metrics(m_Font);
                    qreal bandTopY = ((qreal)h) - metrics.height() - 2 * VER_MARGIN - m_BandPlanHeight;
                    QList<BandInfo> hoverBands = BandPlan::Get().getBandsEncompassing(hoverFrequency);
                    if(m_BandPlanEnabled && py > bandTopY && !hoverBands.empty())
                    {
                        for (auto & hoverBand : hoverBands)
                            toolTipText.append("\n" + hoverBand.name);
                    }
                    showToolTip(event, toolTipText);
                }
            }
            m_GrabPosition = 0;
        }
    }
    else
    {
        // not in Overlay region
        if (event->buttons() == Qt::NoButton)
        {
            if (NOCAP != m_CursorCaptured)
                setCursor(QCursor(Qt::ArrowCursor));

            m_CursorCaptured = NOCAP;
            m_GrabPosition = 0;
        }
        if (m_TooltipsEnabled)
        {
            const quint64 line_ms = msecFromY(py);
            QString timeStr;
            if (line_ms >= wf_valid_since_ms)
            {
                QDateTime tt;
                tt.setMSecsSinceEpoch(msecFromY(py));
                timeStr = tt.toString("yyyy.MM.dd hh:mm:ss.zzz");
            }
            else{
                timeStr = "[time not valid]";
            }

            showToolTip(event, QString("%1\n%2 kHz")
                                       .arg(timeStr)
                                       .arg(freqFromX(px)/1.e3, 0, 'f', 3));
        }
    }
    // process mouse moves while in cursor capture modes
    if (YAXIS == m_CursorCaptured)
    {
        if (event->buttons() & Qt::LeftButton)
        {
            setCursor(QCursor(Qt::ClosedHandCursor));
            // move Y scale up/down
            float delta_px = m_Yzero - py;
            float delta_db = delta_px * fabsf(m_PandMindB - m_PandMaxdB) /
                             (float)h;
            m_PandMindB -= delta_db;
            m_PandMaxdB -= delta_db;
            if (out_of_range(m_PandMindB, m_PandMaxdB))
            {
                m_PandMindB += delta_db;
                m_PandMaxdB += delta_db;
            }
            else
            {
                emit pandapterRangeChanged(m_PandMindB, m_PandMaxdB);

                m_histIIRValid = false;

                m_Yzero = py;

                updateOverlay();
            }
        }
    }
    else if (XAXIS == m_CursorCaptured)
    {
        if (event->buttons() & (Qt::LeftButton | Qt::MiddleButton))
        {
            setCursor(QCursor(Qt::ClosedHandCursor));
            // pan viewable range or move center frequency
            int delta_px = m_Xzero - px;
            qint64 delta_hz = qRound64((qreal)delta_px * (qreal)m_Span / (qreal)w);
            if (delta_hz != 0) // update m_Xzero only on real change
            {
                if (event->buttons() & Qt::MiddleButton)
                {
                    m_CenterFreq += delta_hz;
                    m_DemodCenterFreq += delta_hz;
                    emit newDemodFreq(m_DemodCenterFreq, m_DemodCenterFreq - m_CenterFreq);
                }
                else
                {
                    setFftCenterFreq(m_FftCenter + delta_hz);
                }

                m_MaxHoldValid = false;
                m_MinHoldValid = false;
                m_histIIRValid = false;

                m_Xzero = px;

                updateOverlay();
            }
        }
    }
    else if (LEFT == m_CursorCaptured)
    {
        // moving in demod lowcut region
        if (event->buttons() & (Qt::LeftButton | Qt::RightButton))
        {
            // moving in demod lowcut region with left button held
            if (m_GrabPosition != 0)
            {
                m_DemodLowCutFreq = freqFromX(px - m_GrabPosition ) - m_DemodCenterFreq;
                m_DemodLowCutFreq = std::min(m_DemodLowCutFreq, m_DemodHiCutFreq - FILTER_WIDTH_MIN_HZ);
                m_DemodLowCutFreq = roundFreq(m_DemodLowCutFreq, m_FilterClickResolution);

                if (m_symetric && (event->buttons() & Qt::LeftButton))  // symmetric adjustment
                {
                    m_DemodHiCutFreq = -m_DemodLowCutFreq;
                }
                clampDemodParameters();

                emit newFilterFreq(m_DemodLowCutFreq, m_DemodHiCutFreq);
                updateOverlay();
            }
            else
            {
                // save initial grab position from m_DemodFreqX
                m_GrabPosition = px-m_DemodLowCutFreqX;
            }
        }
        else if (event->buttons() & ~Qt::NoButton)
        {
            setCursor(QCursor(Qt::ArrowCursor));
            m_CursorCaptured = NOCAP;
        }
    }
    else if (RIGHT == m_CursorCaptured)
    {
        // moving in demod highcut region
        if (event->buttons() & (Qt::LeftButton | Qt::RightButton))
        {
            // moving in demod highcut region with right button held
            if (m_GrabPosition != 0)
            {
                m_DemodHiCutFreq = freqFromX(px-m_GrabPosition) - m_DemodCenterFreq;
                m_DemodHiCutFreq = std::max(m_DemodHiCutFreq, m_DemodLowCutFreq + FILTER_WIDTH_MIN_HZ);
                m_DemodHiCutFreq = roundFreq(m_DemodHiCutFreq, m_FilterClickResolution);

                if (m_symetric && (event->buttons() & Qt::LeftButton)) // symmetric adjustment
                {
                    m_DemodLowCutFreq = -m_DemodHiCutFreq;
                }
                clampDemodParameters();

                emit newFilterFreq(m_DemodLowCutFreq, m_DemodHiCutFreq);
                updateOverlay();
            }
            else
            {
                // save initial grab position from m_DemodFreqX
                m_GrabPosition = px - m_DemodHiCutFreqX;
            }
        }
        else if (event->buttons() & ~Qt::NoButton)
        {
            setCursor(QCursor(Qt::ArrowCursor));
            m_CursorCaptured = NOCAP;
        }
    }
    else if (CENTER == m_CursorCaptured)
    {
        // moving inbetween demod lowcut and highcut region
        if (event->buttons() & Qt::LeftButton)
        {   // moving inbetween demod lowcut and highcut region with left button held
            if (m_GrabPosition != 0)
            {
                m_DemodCenterFreq = roundFreq(freqFromX(px - m_GrabPosition),
                                              m_ClickResolution );
                emit newDemodFreq(m_DemodCenterFreq,
                                  m_DemodCenterFreq - m_CenterFreq);
                updateOverlay();
            }
            else
            {
                // save initial grab position from m_DemodFreqX
                m_GrabPosition = px - m_DemodFreqX;
            }
        }
        else if (event->buttons() & ~Qt::NoButton)
        {
            setCursor(QCursor(Qt::ArrowCursor));
            m_CursorCaptured = NOCAP;
        }
    }
    else if (MARKER_A == m_CursorCaptured
             && px < w - m_CursorCaptureDelta
             && px > m_YAxisWidth + m_CursorCaptureDelta)
    {
        if (event->buttons() & Qt::LeftButton)
        {
            qint64 prevA = m_MarkerFreqA;
            m_MarkerFreqA = freqFromX(px);
            emit markerSelectA(m_MarkerFreqA);
            // Shift-drag moves both markers
            if ((event->modifiers() & Qt::ShiftModifier) && m_MarkerFreqB != MARKER_OFF) {
                qint64 df = m_MarkerFreqA - prevA;
                m_MarkerFreqB += df;
                emit markerSelectB(m_MarkerFreqB);
            }
        }
        else if (event->buttons() & ~Qt::NoButton)
        {
            setCursor(QCursor(Qt::ArrowCursor));
            m_CursorCaptured = NOCAP;
        }
    }
    else if (MARKER_B == m_CursorCaptured
             && px < w - m_CursorCaptureDelta
             && px > m_YAxisWidth + m_CursorCaptureDelta)
    {
        if (event->buttons() & Qt::LeftButton)
        {
            qint64 prevB = m_MarkerFreqB;
            m_MarkerFreqB = freqFromX(px);
            emit markerSelectB(m_MarkerFreqB);
            // Shift-drag moves both markers
            if ((event->modifiers() & Qt::ShiftModifier) && m_MarkerFreqA != MARKER_OFF) {
                qint64 df = m_MarkerFreqB - prevB;
                m_MarkerFreqA += df;
                emit markerSelectA(m_MarkerFreqA);
            }
        }
        else if (event->buttons() & ~Qt::NoButton)
        {
            setCursor(QCursor(Qt::ArrowCursor));
            m_CursorCaptured = NOCAP;
        }
    }
    else
    {
        // cursor not captured
        m_GrabPosition = 0;
    }
    if (!this->rect().contains(pt))
    {
        if (NOCAP != m_CursorCaptured)
            setCursor(QCursor(Qt::ArrowCursor));
        m_CursorCaptured = NOCAP;
    }
}


int CPlotter::getNearestPeak(QPoint pt)
{
    int px = qRound((qreal)pt.x() * m_DPR);
    int py = qRound((qreal)pt.y() * m_DPR);

    QMap<int, qreal>::const_iterator i = m_Peaks.lowerBound(px - PEAK_CLICK_MAX_H_DISTANCE);
    QMap<int, qreal>::const_iterator upperBound = m_Peaks.upperBound(px + PEAK_CLICK_MAX_H_DISTANCE);
    qreal   dist = 1.0e10;
    int     best = -1;

    for ( ; i != upperBound; i++)
    {
        int x = i.key();
        qreal y = i.value();

        if (abs(y - py) > PEAK_CLICK_MAX_V_DISTANCE)
            continue;

        qreal d = pow(y - py, 2) + pow(x - px, 2);
        if (d < dist)
        {
            dist = d;
            best = x;
        }
    }

    return best;
}

/** Set waterfall span in milliseconds */
void CPlotter::setWaterfallSpan(quint64 span_ms)
{
    wf_span = span_ms;
    quint64 tnow = QDateTime::currentMSecsSinceEpoch();
    if (!m_WaterfallImage.isNull()) {
        wf_epoch = tnow;
        wf_count = 0;
        msec_per_wfline = (double)wf_span / (qreal)m_WaterfallImage.height();
    }
    wf_valid_since_ms = tnow;
    clearWaterfallBuf();
}

void CPlotter::clearWaterfallBuf()
{
    for (int i = 0; i < MAX_SCREENSIZE; i++)
        m_wfbuf[i] = 0.0;
}

/** Get waterfall time resolution in milleconds / line. */
quint64 CPlotter::getWfTimeRes() const
{
    if (msec_per_wfline)
        return msec_per_wfline;
    else
        // Auto mode, interval is rounded down to nearest int div
        return 1000 / fft_rate;
}

void CPlotter::setFftRate(int rate_hz)
{
    fft_rate = rate_hz;
    wf_valid_since_ms = QDateTime::currentMSecsSinceEpoch();
    clearWaterfallBuf();
}

// Called when a mouse button is pressed
void CPlotter::mousePressEvent(QMouseEvent * event)
{
    QPoint pt = event->pos();
    int px = qRound((qreal)pt.x() * m_DPR);
    int py = qRound((qreal)pt.y() * m_DPR);
    QPoint ppos = QPoint(px, py);

    if (NOCAP == m_CursorCaptured)
    {
        if (isPointCloseTo(px, m_DemodFreqX, m_CursorCaptureDelta))
        {
            // move demod box center frequency region
            m_CursorCaptured = CENTER;
            m_GrabPosition = px - m_DemodFreqX;
        }
        else if (isPointCloseTo(px, m_DemodLowCutFreqX, m_CursorCaptureDelta))
        {
            // filter low cut
            m_CursorCaptured = LEFT;
            m_GrabPosition = px - m_DemodLowCutFreqX;
        }
        else if (isPointCloseTo(px, m_DemodHiCutFreqX, m_CursorCaptureDelta))
        {
            // filter high cut
            m_CursorCaptured = RIGHT;
            m_GrabPosition = px - m_DemodHiCutFreqX;
        }
        else
        {
            if (event->buttons() == Qt::LeftButton)
            {
                // {shift|ctrl|ctrl-shift}-left-click: set ab markers around signal at cursor
                quint32 mods = event->modifiers() & (Qt::ShiftModifier|Qt::ControlModifier);
                if (m_MarkersEnabled && ((event->modifiers() & mods) != 0))
                {
                    float *selectBuf = nullptr;

                    // when max hold is valid, ctrl-shift selects max hold
                    if (m_MaxHoldActive && (mods == (Qt::ShiftModifier | Qt::ControlModifier)))
                    {
                        selectBuf = m_fftMaxHoldBuf;
                    }
                    // in max mode, shift selects max
                    else if (m_PlotMode == PLOT_MODE_MAX && (mods == Qt::ShiftModifier))
                    {
                        selectBuf = m_fftMaxBuf;
                    }
                    // in avg mode, shift select avg
                    else if (m_PlotMode == PLOT_MODE_AVG && (mods == Qt::ShiftModifier))
                    {
                        selectBuf = m_fftAvgBuf;
                    }
                    // in filled and histogram modes, shift selects max, ctrl selects avg
                    else if (m_PlotMode == PLOT_MODE_FILLED || m_PlotMode == PLOT_MODE_HISTOGRAM)
                    {
                        if (mods == Qt::ShiftModifier)
                        {
                            selectBuf = m_fftAvgBuf;
                        }
                        else if (mods == Qt::ControlModifier)
                        {
                            selectBuf = m_fftMaxBuf;
                        }
                    }

                    // ignore if data source is not valid
                    if (m_fftDataSize && selectBuf)
                    {
                        // Find the data value of the click y()

                        const qreal plotHeight = m_2DPixmap.height();
                        const float panddBGainFactor = (float)plotHeight / fabsf(m_PandMaxdB - m_PandMindB);
                        const float vlog = m_PandMaxdB - py / panddBGainFactor;
                        const float v = powf(10.0f, vlog / 10.0f);

                        // Ignore clicks exactly on the plot, below the
                        // pandapter, or when uninitialized
                        if (v != selectBuf[px]
                            && py < plotHeight
                            && m_fftDataSize > 0)
                        {
                            int xLeft = px;
                            int xRight = px;
                            // Select span below the plot
                            if (v < selectBuf[px])
                            {
                                for(; xLeft > 0 && selectBuf[xLeft] > v; --xLeft);
                                for(; xRight < m_fftDataSize && selectBuf[xRight] > v; ++xRight);
                            }
                            // Select span above the plot
                            else
                            {
                                for(; xLeft > 0 && selectBuf[xLeft] < v; --xLeft);
                                for(; xRight < m_fftDataSize && selectBuf[xRight] < v; ++xRight);
                            }
                            qint64 freqLeft = freqFromX(xLeft);
                            qint64 freqRight = freqFromX(xRight);

                            emit markerSelectA(freqLeft);
                            emit markerSelectB(freqRight);
                        }
                    }
                }

                // left-click with no modifiers: set center frequency
                else if (mods == 0) {
                    int best = -1;

                    if (m_PeakDetectActive > 0)
                        best = getNearestPeak(pt);
                    if (best != -1)
                        m_DemodCenterFreq = freqFromX(best);
                    else
                        m_DemodCenterFreq = roundFreq(freqFromX(px), m_ClickResolution);

                    // if cursor not captured set demod frequency and start demod box capture
                    emit newDemodFreq(m_DemodCenterFreq, m_DemodCenterFreq - m_CenterFreq);

                    // save initial grab position from m_DemodFreqX
                    // setCursor(QCursor(Qt::CrossCursor));
                    m_CursorCaptured = CENTER;
                    m_GrabPosition = 1;
                    updateOverlay();
                }
            }
            else if (event->buttons() == Qt::MiddleButton)
            {
                // set center freq
                m_CenterFreq = roundFreq(freqFromX(px), m_ClickResolution);
                m_DemodCenterFreq = m_CenterFreq;
                emit newDemodFreq(m_DemodCenterFreq, m_DemodCenterFreq - m_CenterFreq);
                updateOverlay();
            }
            else if (event->buttons() == Qt::RightButton)
            {
                // reset frequency zoom
                resetHorizontalZoom();
            }
        }
    }
    else
    {
        if (m_CursorCaptured == YAXIS)
            // get ready for moving Y axis
            m_Yzero = py;
        else if (m_CursorCaptured == XAXIS)
        {
            m_Xzero = px;
            if (event->buttons() == Qt::RightButton)
            {
                // reset frequency zoom
                resetHorizontalZoom();
            }
        }
        else if (m_CursorCaptured == TAG)
        {
            for (auto & tag : m_Taglist)
            {
                if (tag.first.contains(ppos))
                {
                    m_DemodCenterFreq = tag.second;
                    emit newDemodFreq(m_DemodCenterFreq, m_DemodCenterFreq - m_CenterFreq);
                    break;
                }
            }
        }
    }
}

void CPlotter::mouseReleaseEvent(QMouseEvent * event)
{
    QPoint pt = event->pos();
    int py = qRound((qreal)pt.y() * m_DPR);

    if (py >= m_OverlayPixmap.height())
    {
        // not in Overlay region
        if (NOCAP != m_CursorCaptured)
            setCursor(QCursor(Qt::ArrowCursor));

        m_CursorCaptured = NOCAP;
        m_GrabPosition = 0;
    }
    else
    {
        if (YAXIS == m_CursorCaptured)
        {
            setCursor(QCursor(Qt::OpenHandCursor));
            m_Yzero = -1;
        }
        else if (XAXIS == m_CursorCaptured)
        {
            setCursor(QCursor(Qt::OpenHandCursor));
            m_Xzero = -1;
        }
    }
}


// Make a single zoom step on the X axis.
void CPlotter::zoomStepX(float step, int x)
{
    // Limit zoom out to 1.0 and zoom in to where there are 5 fft points on the
    // screen. m_fftDataSize is initialized to 0 ... if the app hasn't started
    // yet, allow any zoom level.
    if (m_fftDataSize != 0)
    {
        double currentZoom = (double)m_SampleFreq / (double)m_Span;
        if ((step >= 1.0f && currentZoom <= 1.0)
            || (step < 1.0f && currentZoom >= (double)m_fftDataSize / 4))
            return;
    }

    // calculate new range shown on FFT
    float new_span = std::min(
        (float)m_Span * (float)step, (float)m_SampleFreq);

    // Keep frequency under pointer the same and calculated the offset to the
    // plot center.
    float offset = (float)(freqFromX(x) - m_CenterFreq - m_FftCenter);
    float new_FftCenter = (float)m_FftCenter + offset * (1.0f - step);

    // Keep edges of plot in valid frequency range. The plot may need to be
    // panned.
    const float max_limit = (float)m_SampleFreq / 2.0f;
    const float min_limit = - (float)m_SampleFreq / 2.0f;
    float f_max = new_FftCenter + new_span / 2.0f;
    float f_min = new_FftCenter - new_span / 2.0f;
    if (f_min < min_limit)
    {
        f_min = min_limit;
        f_max = f_min + new_span;
    }
    if (f_max > max_limit)
    {
        f_max = max_limit;
        f_min = f_max - new_span;
    }

    // Make span into an even integer.
    quint32 new_span_int = qRound(f_max - f_min);
    if( new_span_int & 1 )
    {
        new_span_int--;
    }

    // Explicitly set m_Span instead of calling setSpanFreq(), which also calls
    // setFftCenterFreq() and updateOverlay() internally. Span needs to be set
    // before frequency limits can be checked in setFftCenterFreq().
    m_Span = new_span;
    setFftCenterFreq(qRound64((f_max + f_min) / 2.0f));

    m_MaxHoldValid = false;
    m_MinHoldValid = false;
    m_histIIRValid = false;

    updateOverlay();

    double zoom = (double)m_SampleFreq / (double)m_Span;
    emit newZoomLevel(zoom);
    qCDebug(plotter) << QString("Spectrum zoom: %1x").arg(zoom, 0, 'f', 1);
}

// Zoom on X axis (absolute level)
void CPlotter::zoomOnXAxis(float level)
{
    float current_level = (float)m_SampleFreq / (float)m_Span;
    zoomStepX(current_level / level, qRound((qreal)m_Size.width() * m_DPR / 2.0));
    updateOverlay();
}

void CPlotter::setPlotMode(int mode)
{
    m_PlotMode = (ePlotMode)mode;
    m_MaxHoldValid = false;
    m_MinHoldValid = false;
    // Do not need to invalidate IIR data when switching modes

    updateOverlay();
}

void CPlotter::setPlotScale(int scale, bool perHz)
{
    m_PlotScale = (ePlotScale)scale;
    m_PlotPerHz = perHz;
    m_MaxHoldValid = false;
    m_MinHoldValid = false;
    m_IIRValid = false;
    m_histIIRValid = false;
}

void CPlotter::setWaterfallMode(int mode)
{
    m_WaterfallMode = (eWaterfallMode)mode;
}

// Called when a mouse wheel is turned
void CPlotter::wheelEvent(QWheelEvent * event)
{
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
    QPoint pt = QPoint(event->pos());
#else
    QPointF pt = event->position();
#endif
    int h = m_OverlayPixmap.height();
    int px = qRound((qreal)pt.x() * m_DPR);
    int py = qRound((qreal)pt.y() * m_DPR);

    // delta is in eigths of a degree, 15 degrees is one step
    int delta = m_InvertScrolling? -event->angleDelta().y() : event->angleDelta().y();
    double numSteps = delta / (8.0 * 15.0);
    // zoom faster when Ctrl is held
    double zoomBase = (event->modifiers() & Qt::ControlModifier) ? 0.7 : 0.9;

    if (m_CursorCaptured == YAXIS)
    {
        // Vertical zoom. Wheel down: zoom out, wheel up: zoom in
        // During zoom we try to keep the point (dB or kHz) under the cursor fixed
        float zoom_fac = pow(zoomBase, numSteps);
        float ratio = (float) py / (float) h;
        float db_range = m_PandMaxdB - m_PandMindB;
        float y_range = (float) h;
        float db_per_pix = db_range / y_range;
        float fixed_db = m_PandMaxdB - py * db_per_pix;

        db_range = qBound(FFT_MIN_DB_RANGE, db_range * zoom_fac, FFT_MAX_DB - FFT_MIN_DB);
        m_PandMaxdB = fixed_db + ratio * db_range;
        if (m_PandMaxdB > FFT_MAX_DB)
            m_PandMaxdB = FFT_MAX_DB;

        m_PandMindB = m_PandMaxdB - db_range;
        if (m_PandMindB < FFT_MIN_DB)
            m_PandMindB = FFT_MIN_DB;

        m_histIIRValid = false;

        emit pandapterRangeChanged(m_PandMindB, m_PandMaxdB);
    }
    else if (m_CursorCaptured == XAXIS)
    {
        zoomStepX(pow(zoomBase, numSteps), px);
    }
    else if (event->modifiers() & Qt::ControlModifier)
    {
        // filter width
        m_DemodLowCutFreq -= numSteps * m_ClickResolution;
        m_DemodHiCutFreq += numSteps * m_ClickResolution;
        clampDemodParameters();
        emit newFilterFreq(m_DemodLowCutFreq, m_DemodHiCutFreq);
    }

    else if (event->modifiers() & Qt::ShiftModifier)
    {
        // filter shift
        m_DemodLowCutFreq += numSteps * m_ClickResolution;
        m_DemodHiCutFreq += numSteps * m_ClickResolution;
        clampDemodParameters();
        emit newFilterFreq(m_DemodLowCutFreq, m_DemodHiCutFreq);
    }
    else
    {
        // small steps will be lost by roundFreq, let them accumulate
        m_CumWheelDelta += delta;
        if (abs(m_CumWheelDelta) < 8*15)
            return;
        numSteps = m_CumWheelDelta / (8.0 * 15.0);

        // inc/dec demod frequency
        m_DemodCenterFreq += (numSteps * m_ClickResolution);
        m_DemodCenterFreq = roundFreq(m_DemodCenterFreq, m_ClickResolution );
        emit newDemodFreq(m_DemodCenterFreq, m_DemodCenterFreq-m_CenterFreq);
    }

    updateOverlay();
    m_CumWheelDelta = 0;
}

// Called when screen size changes so must recalculate bitmaps
void CPlotter::resizeEvent(QResizeEvent* )
{
    if (!size().isValid())
        return;

    m_DPR = devicePixelRatioF();
    QSize s = QSize(size().width(), size().height());
    if (m_Size != s)
    {
        m_Size = s;

        // Use scaled system font
        m_Font = QFont();
        m_Font.setPointSizeF(m_Font.pointSizeF() * m_DPR);

        // Higher resolution pixmaps are used with higher DPR. They are
        // rescaled in paintEvent().
        const int w = qRound((qreal)s.width() * m_DPR);
        const int plotHeight = qRound((qreal)m_Percent2DScreen * (qreal)s.height() / 100.0 * m_DPR);
        const int wfHeight = qRound((qreal)s.height() * m_DPR) - plotHeight;

        m_OverlayPixmap = QPixmap(w, plotHeight);
        m_OverlayPixmap.fill(Qt::transparent);

        m_2DPixmap = QPixmap(w, plotHeight);
        m_2DPixmap.fill(QColor::fromRgba(PLOTTER_BGD_COLOR));

        // No waterfall, use null image
        if (wfHeight == 0)
        {
            m_WaterfallImage = QImage();
        }

        // New waterfall, create blank area
        else if (m_WaterfallImage.isNull()) {
            m_WaterfallImage = QImage(w, wfHeight, QImage::Format_RGB32);
            m_WaterfallImage.setDevicePixelRatio(m_DPR);
            m_WaterfallImage.fill(Qt::black);
        }

        // Existing waterfall, rescale width but no height as that would
        // invalidate time
        else
        {
            QImage oldWaterfall = m_WaterfallImage.scaled(
                w, m_WaterfallImage.height(),
                Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
            m_WaterfallImage = QImage(w, wfHeight, QImage::Format_RGB32);
            m_WaterfallImage.setDevicePixelRatio(m_DPR);
            m_WaterfallImage.fill(Qt::black);
            memcpy(m_WaterfallImage.bits(), oldWaterfall.bits(),
                m_WaterfallImage.bytesPerLine() * std::min(m_WaterfallImage.height(), oldWaterfall.height()));
        }

        // Invalidate on resize
        m_MaxHoldValid = false;
        m_MinHoldValid = false;
        m_histIIRValid = false;
        // Do not need to invalidate IIR data (just histogram IIR)

        // Waterfall accumulator my be the wrong size now, so invalidate.
        if (msec_per_wfline > 0)
            clearWaterfallBuf();

        // Other things that need to scale with DPR
        m_CursorCaptureDelta = qRound((qreal)CUR_CUT_DELTA * m_DPR);
    }

    updateOverlay();
    emit newSize();
}

// Called by QT when screen needs to be redrawn
void CPlotter::paintEvent(QPaintEvent *)
{
    // Pixmap resolution scales with DPR. Here, they are rescaled to fit the
    // the CPlotter resolution.

    QPainter painter(this);

    int plotHeightT = 0;
    if (!m_2DPixmap.isNull())
    {
        const int plotWidthS = m_2DPixmap.width();
        const int plotHeightS = m_2DPixmap.height();
        const QRectF plotRectS(0.0, 0.0, plotWidthS, plotHeightS);

        const int plotWidthT = qRound((qreal)plotWidthS / m_DPR);
        plotHeightT = qRound((qreal)plotHeightS / m_DPR);
        const QRectF plotRectT(0.0, 0.0, plotWidthT, plotHeightT);

        painter.drawPixmap(plotRectT, m_2DPixmap, plotRectS);
    }

    if (!m_WaterfallImage.isNull())
    {
        painter.drawImage(QPointF(0.0, plotHeightT), m_WaterfallImage);
    }
}

// Called to update spectrum data for displaying on the screen
void CPlotter::draw(bool newData)
{
    qint32        i, j;
    float         histMax;
    QFontMetricsF metrics(m_Font);

    // No fft data yet? Draw overlay if needed and return.
    if (m_fftDataSize == 0)
    {
        if (!m_2DPixmap.isNull()) {
            // Update the overlay if needed
            if (m_DrawOverlay)
            {
                drawOverlay();
                m_DrawOverlay = false;
            }

            // Draw overlay over plot
            m_2DPixmap.fill(QColor::fromRgba(PLOTTER_BGD_COLOR));
            QPainter painter(&m_2DPixmap);
            painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
            painter.drawPixmap(QPointF(0.0, 0.0), m_OverlayPixmap);
            update();
        }

        return;
    }

    QPointF avgLineBuf[MAX_SCREENSIZE];
    QPointF maxLineBuf[MAX_SCREENSIZE];

    const quint64 tnow_ms = QDateTime::currentMSecsSinceEpoch();

    // Pixmaps might be null, so scale up m_Size to get width.
    const qreal w = m_Size.width() * m_DPR;
    const qreal plotHeight = m_2DPixmap.height();
    const qreal shadowOffset = metrics.height() / 20.0;

    // Scale plotter for graph height
    const float panddBGainFactor = (float)plotHeight / fabsf(m_PandMaxdB - m_PandMindB);
    // Scale waterfall and histogram for colormap
    const float wfdBGainFactor = 256.0f / fabsf(m_WfMaxdB - m_WfMindB);

    const double fftSize = m_fftDataSize;
    const double sampleFreq = (double)m_SampleFreq;
    const double fftCenter = (double)m_FftCenter;
    const double span = (double)m_Span;
    const double startFreq = fftCenter - span / 2.0;
    const double binsPerHz = fftSize / sampleFreq;

    // Scale factor for x -> fft bin (pixels per bin).
    double xScale = sampleFreq * w / fftSize / span;

    // Center of fft is the center of the DC bin. The Nyquist bin (index 0
    // after shift) is not used.
    const double startBinD = startFreq * binsPerHz + fftSize / 2.0;
    const qint32 startBin = std::min(qRound(startBinD), m_fftDataSize - 1);
    const qint32 numBins = (qint32)ceil(span * binsPerHz);
    const qint32 endBin = startBin + numBins;
    const qint32 minbin = std::max(startBin, 1);
    const qint32 maxbin = std::min(endBin + 1, m_fftDataSize - 1);

    const qint32 xmin = qRound((double)(minbin - startBin) * xScale);
    const qint32 xmax = std::min(qRound((double)(maxbin - startBin) * xScale), qRound(w));

    const float frameTime = 1.0f / (float)fft_rate;

    // Do plotter work only if visible.
    const bool plotterVisible = (!m_2DPixmap.isNull());

    // Limit plotter drawing rate.
    const bool drawPlotter = (plotterVisible
        && tnow_ms >= tlast_plot_drawn_ms + PLOTTER_UPDATE_LIMIT_MS);

    // Do not waste time with histogram calculations unless in this mode.
    const bool doHistogram = (plotterVisible && m_PlotMode == PLOT_MODE_HISTOGRAM && (!m_histIIRValid || newData));

    // Use fewer histogram bins when statistics are sparse
    const int histBinsDisplayed = std::min(
        MAX_HISTOGRAM_SIZE,
        std::max(32,
            qRound(32 * (float)numBins / 2048.0f))
        );

    // Amount to add to histogram for each hit
    const float histWeight = 10e6f * frameTime / (float)histBinsDisplayed / (float)fftSize;

    // Bins / dB
    const float histdBGainFactor = (float)histBinsDisplayed / fabsf(m_PandMaxdB - m_PandMindB);

    // Show max and average highlights on histogram if it would not be too
    // cluttered
    const bool showHistHighlights = histBinsDisplayed >= MAX_HISTOGRAM_SIZE / 2;

    // Waterfall is advanced only if visible and running, and if there is new
    // data. Repaints for other reasons do not require any action here.
    const bool doWaterfall = !m_WaterfallImage.isNull() && m_Running && newData;

    // Draw avg line, except in max mode. Suppress if it would clutter histogram.
    const bool doAvgLine = m_PlotMode != PLOT_MODE_MAX
                           && (m_PlotMode != PLOT_MODE_HISTOGRAM || showHistHighlights);

    // Draw max line, except in avg and histogram modes
    const bool doMaxLine = m_PlotMode != PLOT_MODE_AVG
                           && m_PlotMode != PLOT_MODE_HISTOGRAM;

    // Initialize results
    if (doHistogram)
        memset(m_histogram, 0, sizeof(m_histogram));

    // Peak means "peak of average" in AVG mode, else "peak of max"
    const bool peakIsAverage = m_PlotMode == PLOT_MODE_AVG;
    // Min mean "min of peak" in PEAK mode, else "min of average"
    const bool minIsAverage = m_PlotMode != PLOT_MODE_MAX;

    // Make sure zeros don't get through to log calcs
    const float fmin = std::numeric_limits<float>::min();

    float vmax;
    float vmaxIIR;
    float vsum;
    float vsumIIR;

    if ((qreal)numBins >= w)
    {
        qint32 count;
        qint32 xprev = xmin;
        bool first = true;

        for(qint32 i = minbin; i <= maxbin; i++)
        {
            const float xD = (float)(i - startBin) * (float)xScale;
            const int x = qRound(xD);

            // Plot uses IIR output. Histogram and waterfall use raw fft data.
            const float v = m_fftData[i];
            const float viir = m_fftIIR[i];

            if (first)
            {
                vmax = v;
                vmaxIIR = viir;
                vsum = v;
                vsumIIR = viir;
                count = 1;
            }

            // Histogram increments the appropriate bin for each value. Ignore
            // out-of-range values, rather than clipping. Allocate value to two
            // closest bins using linear interpolation.
            if (doHistogram)
            {
                const float binD = histdBGainFactor * (m_PandMaxdB - 10.0f * log10f(v));
                if (binD > 0.0f && binD < (float)histBinsDisplayed) {
                    const int binLeft = std::max((int)(xD - 0.5f), 0);
                    const int binRight = std::min(binLeft + 1, numBins - 1);
                    const int binLow = std::min(std::max((int)(binD - 0.5f), 0), histBinsDisplayed - 1);
                    const int binHigh = std::min(binLow + 1, histBinsDisplayed - 1);
                    const float wgtH = (xD - (float)binLeft) / 2.0f;
                    const float wgtV = (binD - (float)binLow) / 2.0f;
                    m_histogram[binLeft][binLow] += (1.0f - wgtV) * (1.0f - wgtH) * histWeight;
                    m_histogram[binLeft][binHigh] += wgtV * (1.0f - wgtH) * histWeight;
                    m_histogram[binRight][binLow] += (1.0f - wgtV) * wgtH * histWeight;
                    m_histogram[binRight][binHigh] += wgtV * wgtH * histWeight;
                }
            }

            // New (or last) pixel - output values
            if (x != xprev || i == maxbin)
            {
                vmax = std::max(vmax, fmin);
                m_wfMaxBuf[xprev] = vmax;

                vmaxIIR = std::max(vmaxIIR, fmin);
                m_fftMaxBuf[xprev] = vmaxIIR;

                const float vavg = std::max((float)(vsum / (float)count), fmin);
                m_wfAvgBuf[xprev] = vavg;
                const float vavgIIR = std::max((float)(vsumIIR / (float)count), fmin);
                m_fftAvgBuf[xprev] = vavgIIR;

                // New peak hold value if greater, or reset
                const float currentPeak = m_fftMaxHoldBuf[xprev];
                const float newPeak = peakIsAverage ? vavgIIR : vmaxIIR;
                m_fftMaxHoldBuf[xprev] = m_MaxHoldValid ? std::max(currentPeak, newPeak) : newPeak;

                // New min hold value if less, or reset
                const float currentMin = m_fftMinHoldBuf[xprev];
                const float newMin = minIsAverage ? vavgIIR : vmaxIIR;
                m_fftMinHoldBuf[xprev] = m_MinHoldValid ? std::min(currentMin, newMin) : newMin;

                vmax = v;
                vmaxIIR = viir;
                vsum = v;
                vsumIIR = viir;
                count = 1;
                xprev = x;
            }

            else if (!first)
            {
                vmax = std::max(v, vmax);
                vmaxIIR = std::max(viir, vmaxIIR);
                vsum += v;
                vsumIIR += viir;
                ++count;
            }

            first = false;
        }

        m_MaxHoldValid = true;
        m_MinHoldValid = true;
    }
    // w > m_fftDataSize uses no averaging
    else
    {
        for (i = xmin; i < xmax; i++)
        {
            j = qRound((float)i / (float)xScale + (float)startBinD);

            const float v = m_fftData[j];
            const float viir = m_fftIIR[j];

            m_wfMaxBuf[i] = v;
            m_wfAvgBuf[i] = v;
            m_fftMaxBuf[i] = viir;
            m_fftAvgBuf[i] = viir;

            // New peak hold value if greater, or reset
            const float currentPeak = m_fftMaxHoldBuf[i];
            m_fftMaxHoldBuf[i] = m_MaxHoldValid ? std::max(currentPeak, viir) : viir;

            // New min hold value if less, or reset
            const float currentMin = m_fftMinHoldBuf[i];
            m_fftMinHoldBuf[i] = m_MinHoldValid ? std::min(currentMin, viir) : viir;

            // Histogram increments the appropriate bin for each value. Ignore
            // out-of-range values, rather than clipping. Allocate value to two
            // closest bins using linear interpolation.
            if (doHistogram)
            {
                const float binD = histdBGainFactor * (m_PandMaxdB - 10.0f * log10f(v));
                if (binD > 0.0f && binD < (float)histBinsDisplayed) {
                    const int binLow = std::min(std::max((int)(binD - 0.5f), 0), histBinsDisplayed - 1);
                    const int binHigh = std::min(binLow + 1, histBinsDisplayed - 1);
                    const float wgt = (binD - (float)binLow) / 2.0f;
                    m_histogram[i][binLow] += (1.0f - wgt) * histWeight;
                    m_histogram[i][binHigh] += wgt * histWeight;
                }
            }
        }
    }

    const int npts = xmax - xmin;

    if (doWaterfall)
    {
        // Pick max or avg for waterfall
        float *dataSource;
        if (m_WaterfallMode == WATERFALL_MODE_AVG)
        {
            dataSource = m_wfAvgBuf;
        }
        else if (m_WaterfallMode == WATERFALL_MODE_SYNC)
        {
            if (m_PlotMode == PLOT_MODE_MAX)
            {
                dataSource = m_fftMaxBuf;
            }
            else {
                dataSource = m_fftAvgBuf;
            }
        }
        else // WATERFALL_MODE_MAX
        {
            dataSource = m_wfMaxBuf;
        }

        // if not in "auto" mode, store max waterfall data in accumulator
        if (msec_per_wfline > 0)
        {
            // In avg mode, accumulate so average of frames can be shown
            if (m_WaterfallMode != WATERFALL_MODE_MAX)
            {
                ++wf_avg_count;
                for (i = 0; i < npts; ++i)
                    m_wfbuf[i] += dataSource[i];
            }
            // In max mode, track the max bin over time
            else
            {
                for (i = 0; i < npts; ++i)
                    m_wfbuf[i] = std::max(m_wfbuf[i], dataSource[i]);
            }
        }

        // is it time to update waterfall? msec_per_wfline is 0 in auto mode.
        if (tnow_ms - wf_epoch > wf_count * msec_per_wfline)
        {
            ++wf_count;

            // cursor times are relative to last time drawn
            tlast_wf_ms = tnow_ms;
            if (wf_valid_since_ms == 0)
                wf_valid_since_ms = tnow_ms;
            tlast_wf_drawn_ms = tnow_ms;

            // move current data down one line(must do before attaching a QPainter object)
            memmove(m_WaterfallImage.scanLine(1), m_WaterfallImage.scanLine(0),
                m_WaterfallImage.bytesPerLine() * (m_WaterfallImage.height() - 1));

            // draw new line of fft data at top of waterfall bitmap
            // draw black areas where data will not be draw
            memset(m_WaterfallImage.scanLine(0), 0, m_WaterfallImage.bytesPerLine());

            const bool useWfBuf = msec_per_wfline > 0;
            float _lineFactor;
            if (useWfBuf && m_WaterfallMode != WATERFALL_MODE_MAX)
                _lineFactor = 1.0f / (float)wf_avg_count;
            else
                _lineFactor = 1.0f;
            const float lineFactor = _lineFactor;
            wf_avg_count = 0;

            // Use buffer (max or average) if in manual mode, else current data
            for (i = 0; i < npts; ++i)
            {
                const int ix = i + xmin;
                const float v = useWfBuf ? m_wfbuf[ix] * lineFactor : dataSource[ix];
                qint32 cidx = qRound((m_WfMaxdB - 10.0f * log10f(v)) * wfdBGainFactor);
                cidx = std::max(std::min(cidx, 255), 0);
                m_WaterfallImage.setPixel(ix, 0, m_ColorTbl[255 - cidx].rgb());
            }

            wf_avg_count = 0;
            if (msec_per_wfline > 0)
                clearWaterfallBuf();
        }
    }

    // Update histogram IIR if it will be used.
    if (doHistogram)
    {
        const float gamma = 1.0f;
        const float a = powf(1.0f - m_alpha, gamma);
        // fast attack ... leaving alternative here in case it's useful
        const float aAttack = 1.0f;
        // const float aAttack = 1.0 - a * frameTime;
        const float aDecay = 1.0f - powf(a, 4.0f * frameTime);

        histMax = 0.0;
        for (i = xmin; i < xmax; ++i) {
            for (j = 0; j < histBinsDisplayed; ++j)
            {
                float histV;
                const float histPrev = m_histIIR[i][j];
                const float histNew = m_histogram[i][j];
                // Fast response when invalid
                if (!m_histIIRValid)
                    histV = histNew;
                else
                    histV = histPrev + aAttack * histNew - aDecay * histPrev;
                m_histIIR[i][j] = std::max(histV, 0.0f);
                histMax = std::max(histMax, histV);
            }
        }
        m_histIIRValid = true;

        // 5 Hz time constant for colormap adjustment
        const float histMaxAlpha = std::min(5.0f * frameTime, 1.0f);
        m_histMaxIIR = m_histMaxIIR * (1.0f - histMaxAlpha) + histMax * histMaxAlpha;
    }

    // get/draw the 2D spectrum
    if (drawPlotter)
    {
        tlast_plot_drawn_ms = tnow_ms;

        m_2DPixmap.fill(QColor::fromRgba(PLOTTER_BGD_COLOR));
        QPainter painter2(&m_2DPixmap);
        painter2.translate(QPointF(0.5, 0.5));


        // draw the pandapter
        QBrush fillBrush = QBrush(m_FftFillCol);

        // Fill between max and avg
        QColor maxFillCol = m_FftFillCol;
        maxFillCol.setAlpha(80);
        QBrush maxFillBrush = QBrush(maxFillCol);

        // Diagonal fill for area between markers. Scale the pattern to DPR.
        QColor abFillColor = QColor::fromRgba(PLOTTER_MARKER_COLOR);
        abFillColor.setAlpha(128);
        QBrush abFillBrush = QBrush(abFillColor, Qt::BDiagPattern);

        QColor maxLineColor = QColor(m_FftFillCol);
        if (m_PlotMode == PLOT_MODE_FILLED)
            maxLineColor.setAlpha(128);
        else
            maxLineColor.setAlpha(255);

        QPen maxLinePen = QPen(maxLineColor);

        // Same color as max in avg mode, different for filled mode
        QPen avgLinePen;
        if (m_PlotMode == PLOT_MODE_AVG || m_PlotMode == PLOT_MODE_HISTOGRAM)
        {
            QColor avgLineCol = m_FftFillCol;
            avgLineCol.setAlpha(255);
            avgLinePen = QPen(avgLineCol);
        }
        else {
            QColor avgLineCol = QColor(Qt::cyan);
            avgLineCol.setAlpha(192);
            avgLinePen = QPen(avgLineCol);
        }

        // The m_Marker{AB}X values are one cycle old, which makes for a laggy
        // effect, so get fresh values here.
        const int ax = xFromFreq(m_MarkerFreqA);
        const int bx = xFromFreq(m_MarkerFreqB);
        bool fillMarkers = (m_MarkersEnabled && m_MarkerFreqA != MARKER_OFF
                                             && m_MarkerFreqB != MARKER_OFF);
        const int minMarker = std::min(ax, bx);
        const int maxMarker = std::max(ax, bx);

        const float binSizeY = (float)plotHeight / (float)histBinsDisplayed;
        for (i = 0; i < npts; i++)
        {
            const int ix = i + xmin;
            const qreal ixPlot = (qreal)ix;
            const qreal yMaxD = (qreal)std::max(std::min(
                panddBGainFactor * (m_PandMaxdB - 10.0f * log10f(m_fftMaxBuf[ix])),
                (float)plotHeight), 0.0f);
            const qreal yAvgD = (qreal)std::max(std::min(
                panddBGainFactor * (m_PandMaxdB - 10.0f * log10f(m_fftAvgBuf[ix])),
                (float)plotHeight), 0.0f);

            if (m_PlotMode == PLOT_MODE_HISTOGRAM)
            {
                const float *histData = m_histIIR[(ix)];
                qreal topBin = plotHeight;
                for (j = 0; j < histBinsDisplayed; ++j)
                {
                    qint32 cidx = qRound(histData[j] / m_histMaxIIR * 255.0f * .7f);
                    if (cidx > 0) {
                        cidx += 65;  // 255 * 0.7 = 178, + 65 = 243
                        // Histogram IIR can cause out-of-range cidx
                        cidx = std::max(std::min(cidx, 255), 0);
                        QColor c = m_ColorTbl[cidx];
                        // Paint rectangle
                        const qreal binY = (qreal)binSizeY * j;
                        topBin = std::min(topBin, binY);
                        const qreal binH = (qreal)binSizeY * (j + 1) - binY;
                        painter2.fillRect(QRectF(ixPlot, binY, 1.0, binH), c);
                    }
                }
                // Highlight the top bin, if it isn't too crowded
                if (topBin != plotHeight && showHistHighlights) {
                    painter2.fillRect(QRectF(ixPlot, topBin, 1.0, (qreal)binSizeY), maxLineColor);
                }
            }

            // Add max, average points if they will be drawn
            if (doMaxLine)
                maxLineBuf[i] = QPointF(ixPlot, yMaxD);
            if (doAvgLine)
                avgLineBuf[i] = QPointF(ixPlot, yAvgD);

            // Fill area between markers, even if they are off screen
            qreal yFill = m_PlotMode == PLOT_MODE_MAX ? yMaxD : yAvgD;
            if (fillMarkers && (ix) > minMarker && (ix) < maxMarker) {
                painter2.fillRect(QRectF(ixPlot, yFill + 1.0, 1.0, plotHeight - yFill), abFillBrush);
            }
            if (m_FftFill && m_PlotMode != PLOT_MODE_HISTOGRAM)
            {
                painter2.fillRect(QRectF(ixPlot, yFill + 1.0, 1.0, plotHeight - yFill), m_FftFillCol);
            }
            if (m_PlotMode == PLOT_MODE_FILLED)
            {
                painter2.fillRect(QRectF(ixPlot, yMaxD + 1.0, 1.0, yAvgD - yMaxD), maxFillBrush);
            }
        }

        if (doMaxLine) {
            // NOT scaling to DPR due to performance
            painter2.setPen(maxLinePen);
            painter2.drawPolyline(maxLineBuf, npts);
        }
        if (doAvgLine) {
            // NOT scaling to DPR due to performance
            painter2.setPen(avgLinePen);
            painter2.drawPolyline(avgLineBuf, npts);
        }

        // Max hold
        if (m_MaxHoldActive)
        {
            // Show max(max) except when showing only avg on screen
            for (i = 0; i < npts; i++)
            {
                const int ix = i + xmin;
                const qreal ixPlot = (qreal)ix;
                const qreal yMaxHoldD = (qreal)std::max(std::min(
                    panddBGainFactor * (m_PandMaxdB - 10.0f * log10f(m_fftMaxHoldBuf[ix])),
                    (float)plotHeight), 0.0f);
                maxLineBuf[i] = QPointF(ixPlot, yMaxHoldD);
            }
            // NOT scaling to DPR due to performance
            painter2.setPen(m_MaxHoldColor);
            painter2.drawPolyline(maxLineBuf, npts);

            m_MaxHoldValid = true;
        }

        // Min hold
        if (m_MinHoldActive)
        {
            // Show min(avg) except when showing only max on scree
            for (i = 0; i < npts; i++)
            {
                const int ix = i + xmin;
                const qreal ixPlot = (qreal)ix;
                const qreal yMinHoldD = (qreal)std::max(std::min(
                    panddBGainFactor * (m_PandMaxdB - 10.0f * log10f(m_fftMinHoldBuf[ix])),
                    (float)plotHeight), 0.0f);
                maxLineBuf[i] = QPointF(ixPlot, yMinHoldD);
            }
            // NOT scaling to DPR due to performance
            painter2.setPen(m_MinHoldColor);
            painter2.drawPolyline(maxLineBuf, npts);

            m_MinHoldValid = true;
        }

        // Peak detection
        if (m_PeakDetectActive)
        {
            const int pw = PEAK_WINDOW_HALF_WIDTH;

            // Use data source appropriate for current display mode
            float *_detectSource;
            if (m_MaxHoldActive)
                _detectSource = m_fftMaxHoldBuf;
            else if (m_PlotMode == PLOT_MODE_AVG)
                _detectSource = m_fftAvgBuf;
            else
                _detectSource = m_fftMaxBuf;
            const float *detectSource = _detectSource;

            // Run peak detection periodically. If overlay will be redrawn, run
            // peak detection since zoom/pan may have changed.
            if (tnow_ms > tlast_peaks_ms + PEAK_UPDATE_PERIOD || m_DrawOverlay) {
                tlast_peaks_ms = tnow_ms;
                m_Peaks.clear();

                // Narrow peaks
                for (i = pw; i < npts - pw; ++i) {
                    const int ix = i + xmin;
                    const float vi = detectSource[ix];
                    float sumV = 0;
                    float minV = vi;
                    float maxV = 0;
                    for (j = -pw; j <= pw; ++j) {
                        const float vj = detectSource[ix + j];
                        minV = std::min(minV, vj);
                        maxV = std::max(maxV, vj);
                        sumV += vj;
                    }
                    const float avgV = sumV / (float)(pw * 2 + 1);
                    m_peakSmoothBuf[ix] = avgV;
                    if (vi == maxV && (vi > 2.0f * avgV) && (vi > 4.0f * minV))
                    {
                        const qreal y = (qreal)std::max(std::min(
                            panddBGainFactor * (m_PandMaxdB - 10.0f * log10f(vi)),
                            (float)plotHeight - 0.0f), 0.0f);
                        m_Peaks[ix] = y;
                    }
                }

                // Use the smoothed curve to find wider peaks
                const int pw2 = pw * 5;
                for (i = pw2; i < npts - pw2; ++i) {
                    const int ix = i + xmin;
                    const float vi = m_peakSmoothBuf[ix];
                    float sumV = 0;
                    float minV = vi;
                    float maxV = 0;
                    for (j = -pw2; j <= pw2; ++j) {
                        const float vj = m_peakSmoothBuf[ix + j];
                        minV = std::min(minV, vj);
                        maxV = std::max(maxV, vj);
                        sumV += vj;
                    }
                    const float avgV = sumV / (float)(pw2 * 2);
                    if (vi == maxV && (vi > 2.0f * avgV) && (vi > 4.0f * minV))
                    {
                        const qreal y = (qreal)std::max(std::min(
                            panddBGainFactor * (m_PandMaxdB - 10.0f * log10f(vi)),
                            (float)plotHeight - 0.0f), 0.0f);

                        // Show the wider peak only if there is no very close narrow peak
                        bool found = false;
                        for (j = -pw; j <= pw; ++j) {
                            auto it = m_Peaks.find(ix + j);
                            if (it != m_Peaks.end()) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            m_Peaks[ix] = y;
                        }
                    }
                }
            }

            // Paint peaks with shadow
            QPen peakPen(m_maxFftColor, m_DPR);
            QPen peakShadowPen(Qt::black, m_DPR);
            peakPen.setWidthF(m_DPR);
            for(auto peakx : m_Peaks.keys()) {
                const qreal peakxPlot = (qreal)peakx;
                const qreal peakv = m_Peaks.value(peakx);
                painter2.setPen(peakShadowPen);
                painter2.drawEllipse(
                    QRectF(peakxPlot - 5.0 * m_DPR + shadowOffset,
                           peakv - 5.0 * m_DPR + shadowOffset,
                           10.0 * m_DPR, 10.0 * m_DPR));
                painter2.setPen(peakPen);
                painter2.drawEllipse(
                    QRectF(peakxPlot - 5.0 * m_DPR,
                           peakv - 5.0 * m_DPR,
                           10.0 * m_DPR, 10.0 * m_DPR));
            }
        }

        // Update the overlay if needed
        if (m_DrawOverlay)
        {
            drawOverlay();
            m_DrawOverlay = false;
        }

        // Draw overlay over plot
        painter2.setCompositionMode(QPainter::CompositionMode_SourceOver);
        painter2.drawPixmap(QPointF(0.0, 0.0), m_OverlayPixmap);
    }

    // trigger a new paintEvent
    update();
}

void CPlotter::setRunningState(bool running)
{
    // Reset waterfall time and clear waterfall, since time is no longer correct
    if (running && !m_Running)
    {
        setWaterfallSpan(wf_span);

        // Invalidate any existing data
        m_MaxHoldValid = false;
        m_MinHoldValid = false;
        m_IIRValid = false;
        m_histIIRValid = false;
        m_histMaxIIR = std::numeric_limits<float>::min();

    }

    m_Running = running;
}

/**
 * Set new FFT data.
 * @param fftData Pointer to the new FFT data (same data for pandapter and waterfall).
 * @param size The FFT size.
 *
 * When FFT data is set using this method, the same data will be used for both the
 * pandapter and the waterfall.
 */
void CPlotter::setNewFftData(const float *fftData, int size)
{
    // Make sure zeros don't get through to log calcs
    const float fmin = 1e-20;

    if (size != m_fftDataSize)
    {
        // Reallocate and invalidate IIRs
        m_fftData.resize(size);
        m_fftIIR.resize(size);
        m_X.resize(size);

        m_MaxHoldValid = false;
        m_MinHoldValid = false;
        m_IIRValid = false;

        m_histIIRValid = false;
        m_histMaxIIR = std::numeric_limits<float>::min();

        m_fftDataSize = size;

        // Zoom out if needed to keep about 4 points on the screen
        double currentZoom = (double)m_SampleFreq / (double)m_Span;
        double maxZoom = (double)m_fftDataSize / 4.0;
        if (currentZoom > maxZoom)
            zoomStepX(currentZoom / maxZoom, qRound((qreal)m_Size.width() * m_DPR / 2.0));
    }

    // For dBFS, define full scale as peak (not RMS). A 1.0 FS peak sine wave
    // is 0 dBFS.
    float _pwr_scale = 1.0f / ((float)size * (float)size);

    // For V, convert peak to RMS (/2). 1V peak corresponds to -3.01 dBV (RMS
    // value is 0.707 * peak).
    if (m_PlotScale == PLOT_SCALE_DBV)
        _pwr_scale *= 1.0f / 2.0f;

    // For dBm, the scale is interpreted as V. A 1V peak sine corresponds to
    // 10mW, or 10 dBm. The factor of 2 converts Vpeak to Vrms.
    else if (m_PlotScale == PLOT_SCALE_DBMW50)
        _pwr_scale *= 1000.0f / (2.0f * 50.0f);

    // For units of /Hz, rescale by 1/RBW. For V, this results in /sqrt(Hz), and is
    // used for noise spectral density.
    if (m_PlotPerHz && m_PlotScale != PLOT_SCALE_DBFS)
        _pwr_scale *= (float)size / (float)m_SampleFreq;

    const float pwr_scale = _pwr_scale;
    for (int i = 0; i < size; ++i)
        m_fftData[i] = std::max(fftData[i] * pwr_scale, fmin);

    // Update IIR. If IIR is invalid, set alpha to use latest value. Since the
    // IIR is linear data and users would like to see symmetric attack/decay on
    // the logarithmic y-axis, IIR is in terms of multiplication rather than
    // addition.

    // Time constant, taking update rate into account. Attack and decay rate of
    // change in dB/sec should not visibly change with FFT rate.
    const float _a = powf((float)fft_rate, -1.75f * (1.0f - m_alpha));

    // Make the slider vs alpha nonlinear
    const float gamma = 0.7;
    const float a = powf(_a, gamma);

    // Shortcut expensive pow() if not needed
    const bool needIIR = m_IIRValid                         // Initializing
                      && a != 1.0f;                         // IIR is NOP

    if (needIIR) {
        volk_32f_x2_divide_32f(m_X.data(), m_fftData.data(), m_fftIIR.data(), size);
        volk_32f_s32f_power_32f(m_X.data(), m_X.data(), a, size);
        volk_32f_x2_multiply_32f(m_fftIIR.data(), m_fftIIR.data(), m_X.data(), size);
    }
    else
    {
        memcpy(m_fftIIR.data(), m_fftData.data(), size * sizeof(float));
    }

    m_IIRValid = true;

    draw(true);
}

void CPlotter::setFftAvg(float avg)
{
    m_alpha = avg;
}

void CPlotter::setFftRange(float min, float max)
{
    setWaterfallRange(min, max);
    setPandapterRange(min, max);
}

void CPlotter::setPandapterRange(float min, float max)
{
    if (out_of_range(min, max))
        return;

    m_PandMindB = min;
    m_PandMaxdB = max;
    m_histIIRValid = false;
    updateOverlay();
}

void CPlotter::setWaterfallRange(float min, float max)
{
    if (out_of_range(min, max))
        return;

    m_WfMindB = min;
    m_WfMaxdB = max;
    // no overlay change is necessary
}

// Called to draw an overlay bitmap containing grid and text that
// does not need to be recreated every fft data update.
void CPlotter::drawOverlay()
{
    if (m_OverlayPixmap.isNull())
        return;

    int     x;
    qreal   pixperdiv;
    qreal   adjoffset;
    qreal   dbstepsize;
    qreal   mindbadj;
    QFontMetricsF metrics(m_Font);
    const qreal shadowOffset = metrics.height() / 20.0;
    qreal   w = m_OverlayPixmap.width();
    qreal   h = m_OverlayPixmap.height();

    m_OverlayPixmap.fill(Qt::transparent);
    QPainter painter(&m_OverlayPixmap);
    painter.translate(QPointF(-0.5, -0.5));
    // painter.setRenderHint(QPainter::Antialiasing);
    painter.setFont(m_Font);

    QList<BookmarkInfo> tags;

    // X and Y axis areas
    m_YAxisWidth = metrics.boundingRect("-120").width() + 2 * HOR_MARGIN;
    m_XAxisYCenter = h - metrics.height()/2;
    qreal xAxisHeight = metrics.height() + 2 * VER_MARGIN;
    qreal xAxisTop = h - xAxisHeight;
    qreal fLabelTop = xAxisTop + VER_MARGIN;

    if (m_BookmarksEnabled || m_DXCSpotsEnabled)
    {
        m_Taglist.clear();
        static const QFontMetricsF fm(painter.font());
        static const qreal fontHeight = fm.ascent() + 1;
        static const qreal slant = 5;
        static const qreal levelHeight = fontHeight + 5;
        static const qreal nLevels = h / (levelHeight + slant);
        if (m_BookmarksEnabled)
        {
            tags = Bookmarks::Get().getBookmarksInRange(m_CenterFreq + m_FftCenter - m_Span / 2,
                                                        m_CenterFreq + m_FftCenter + m_Span / 2);
        }
        else
        {
            tags.clear();
        }
        if (m_DXCSpotsEnabled)
        {
            QList<DXCSpotInfo> dxcspots = DXCSpots::Get().getDXCSpotsInRange(m_CenterFreq + m_FftCenter - m_Span / 2,
                                                                             m_CenterFreq + m_FftCenter + m_Span / 2);
            QListIterator<DXCSpotInfo> iter(dxcspots);
            while(iter.hasNext())
            {
                BookmarkInfo tempDXCSpot;
                DXCSpotInfo IterDXCSpot = iter.next();
                tempDXCSpot.name = IterDXCSpot.name;
                tempDXCSpot.frequency = IterDXCSpot.frequency;
                tags.append(tempDXCSpot);
            }
            std::stable_sort(tags.begin(),tags.end());
        }
        QVector<int> tagEnd(nLevels + 1);
        for (auto & tag : tags)
        {
            x = xFromFreq(tag.frequency);
            qreal nameWidth = fm.boundingRect(tag.name).width();

            int level = 0;
            while(level < nLevels && tagEnd[level] > x)
                level++;

            if(level >= nLevels)
            {
                level = 0;
                if (tagEnd[level] > x)
                    continue; // no overwrite at level 0
            }

            tagEnd[level] = x + nameWidth + slant - 1;

            const auto levelNHeight = level * levelHeight;
            const auto levelNHeightBottom = levelNHeight + fontHeight;
            const auto levelNHeightBottomSlant = levelNHeightBottom + slant;

            m_Taglist.append(qMakePair(QRectF(x, levelNHeight, nameWidth + slant, fontHeight), tag.frequency));

            QColor color = QColor(tag.GetColor());
            color.setAlpha(100);
            // Vertical line
            painter.setPen(QPen(color, m_DPR, Qt::DashLine));
            painter.drawLine(QPointF(x, levelNHeightBottomSlant), QPointF(x, xAxisTop));

            // Horizontal line
            painter.setPen(QPen(color, m_DPR, Qt::SolidLine));
            painter.drawLine(QPointF(x + slant, levelNHeightBottom),
                             QPointF(x + nameWidth + slant - 1,
                             levelNHeightBottom));
            // Diagonal line
            painter.drawLine(QPointF(x + 1, levelNHeightBottomSlant - 1),
                             QPointF(x + slant - 1, levelNHeightBottom + 1));

            color.setAlpha(255);
            painter.setPen(QPen(color, 2.0 * m_DPR, Qt::SolidLine));
            painter.drawText(x + slant, levelNHeight, nameWidth,
                             fontHeight, Qt::AlignVCenter | Qt::AlignHCenter,
                             tag.name);
        }
    }

    if (m_BandPlanEnabled)
    {
        QList<BandInfo> bands = BandPlan::Get().getBandsInRange(m_CenterFreq + m_FftCenter - m_Span / 2,
                                                                m_CenterFreq + m_FftCenter + m_Span / 2);

        m_BandPlanHeight = metrics.height() + VER_MARGIN;
        for (auto & band : bands)
        {
            int band_left = std::max(xFromFreq(band.minFrequency), 0);
            int band_right = std::min(xFromFreq(band.maxFrequency), (int)w);
            int band_width = band_right - band_left;
            QRectF rect(band_left, xAxisTop - m_BandPlanHeight, band_width, m_BandPlanHeight);
            painter.fillRect(rect, band.color);
            QString band_label = metrics.elidedText(band.name + " (" + band.modulation + ")", Qt::ElideRight, band_width - 10);
            QRectF textRect(band_left, xAxisTop - m_BandPlanHeight, band_width, metrics.height());
            painter.setPen(QPen(QColor::fromRgba(PLOTTER_TEXT_COLOR), m_DPR));
            painter.drawText(textRect, Qt::AlignCenter, band_label);
        }
    }

    if (m_CenterLineEnabled)
    {
        x = xFromFreq(m_CenterFreq);
        painter.setPen(QPen(QColor::fromRgba(PLOTTER_CENTER_LINE_COLOR), m_DPR));
        painter.drawLine(QPointF(x, 0), QPointF(x, xAxisTop));
    }

    if (m_MarkersEnabled)
    {
        QBrush brush;
        brush.setColor(QColor::fromRgba(PLOTTER_MARKER_COLOR));
        brush.setStyle(Qt::SolidPattern);
        painter.setPen(QPen(QColor::fromRgba(PLOTTER_MARKER_COLOR), m_DPR));

        qreal markerSize = metrics.height() / 2;

        if (m_MarkerFreqA != MARKER_OFF) {
            x = xFromFreq(m_MarkerFreqA);
            m_MarkerAX = x;
            QPolygon poly;
            QPainterPath path;
            poly << QPoint(x - markerSize/2, 0)
                    << QPoint(x + markerSize/2, 0)
                    << QPoint(x, markerSize);
            path.addPolygon(poly);
            painter.drawPolygon(poly);
            painter.fillPath(path, brush);
            painter.drawLine(x, markerSize, x, xAxisTop);
            painter.drawStaticText(QPointF(x + markerSize/2, 0), QStaticText("A"));
        }

        if (m_MarkerFreqB != MARKER_OFF) {
            x = xFromFreq(m_MarkerFreqB);
            m_MarkerBX = x;
            QPolygon poly;
            QPainterPath path;
            poly << QPoint(x - markerSize/2, 0)
                    << QPoint(x + markerSize/2, 0)
                    << QPoint(x, markerSize);
            path.addPolygon(poly);
            painter.drawPolygon(poly);
            painter.fillPath(path, brush);
            painter.drawLine(x, markerSize, x, xAxisTop);
            painter.drawStaticText(QPointF(x + markerSize/2, 0), QStaticText("B"));
        }
    }

    // Frequency grid
    qint64  StartFreq = m_CenterFreq + m_FftCenter - m_Span / 2;
    QString label;
    label.setNum(float((StartFreq + m_Span) / m_FreqUnits), 'f', m_FreqDigits);
    calcDivSize(StartFreq, StartFreq + m_Span,
                qMin(w / (metrics.boundingRect(label).width() + metrics.boundingRect("O").width()),
                     (qreal)HORZ_DIVS_MAX),
                m_StartFreqAdj, m_FreqPerDiv, m_HorDivs);
    pixperdiv = w * (qreal) m_FreqPerDiv / (qreal) m_Span;
    adjoffset = pixperdiv * (qreal) (m_StartFreqAdj - StartFreq) / (qreal) m_FreqPerDiv;

    // Hairline for grid lines
    painter.setPen(QPen(QColor::fromRgba(PLOTTER_GRID_COLOR), 0.0, Qt::DotLine));
    for (int i = 0; i <= m_HorDivs; i++)
    {
        qreal xD = (double)i * pixperdiv + adjoffset;
        if (xD > m_YAxisWidth)
            painter.drawLine(xD, 0, xD, xAxisTop);
    }

    // draw frequency values (x axis)
    makeFrequencyStrs();
    for (int i = 0; i <= m_HorDivs; i++)
    {
        qreal xD = (qreal)i * pixperdiv + adjoffset;
        if (xD > m_YAxisWidth)
        {
            // Shadow
            QRectF shadowRect(xD + shadowOffset - w/2, fLabelTop + shadowOffset,
                              w, metrics.height());
            painter.setPen(QPen(QColor(Qt::black)));
            painter.drawText(shadowRect, Qt::AlignHCenter|Qt::AlignBottom, m_HDivText[i]);
            // Foreground
            QRectF textRect(xD - w/2, fLabelTop,
                            w, metrics.height());
            painter.setPen(QPen(QColor::fromRgba(PLOTTER_TEXT_COLOR)));
            painter.drawText(textRect, Qt::AlignHCenter|Qt::AlignBottom, m_HDivText[i]);
        }
    }

    // Level grid
    qint64 mindBAdj64 = 0;
    qint64 dbDivSize = 0;
    qint64 dbSpan = (qint64) (m_PandMaxdB - m_PandMindB);

    calcDivSize((qint64) m_PandMindB, ((qint64) m_PandMindB) + dbSpan,
                qMax(h / (m_VdivDelta * m_DPR), (qreal)VERT_DIVS_MIN),
                mindBAdj64, dbDivSize, m_VerDivs);

    dbstepsize = (qreal) dbDivSize;
    mindbadj = mindBAdj64;

    pixperdiv = h * (qreal) dbstepsize / (qreal) (m_PandMaxdB - m_PandMindB);
    adjoffset = h * (mindbadj - (qreal) m_PandMindB) / (qreal) (m_PandMaxdB - m_PandMindB);

    qCDebug(plotter) << "minDb =" << m_PandMindB << "maxDb =" << m_PandMaxdB
                     << "mindbadj =" << mindbadj << "dbstepsize =" << dbstepsize
                     << "pixperdiv =" << pixperdiv << "adjoffset =" << adjoffset;

    // Hairline for grid lines
    painter.setPen(QPen(QColor::fromRgba(PLOTTER_GRID_COLOR), 0.0, Qt::DotLine));
    for (int i = 0; i <= m_VerDivs; i++)
    {
        qreal y = h - ((double)i * pixperdiv + adjoffset);
        if (y < h - xAxisHeight)
            painter.drawLine(m_YAxisWidth, y, w, y);
    }

    // draw amplitude values (y axis)
    for (int i = 0; i <= m_VerDivs; i++)
    {
        qreal y = h - ((double)i * pixperdiv + adjoffset);
        qreal th = metrics.height();
        qreal shadowOffset = th / 20.0;
        if ((y < h - xAxisHeight) && (y > th / 2))
        {
            int dB = mindbadj + dbstepsize * i;
            // Shadow
            painter.setPen(QPen(QColor(Qt::black)));
            QRectF shadowRect(HOR_MARGIN + shadowOffset, y - th / 2 + shadowOffset,
                              m_YAxisWidth - 2 * HOR_MARGIN, th);
            painter.drawText(shadowRect, Qt::AlignRight|Qt::AlignVCenter, QString::number(dB));
            // Foreground
            painter.setPen(QPen(QColor::fromRgba(PLOTTER_TEXT_COLOR)));
            QRectF textRect(HOR_MARGIN, y - th / 2,
                            m_YAxisWidth - 2 * HOR_MARGIN, th);
            painter.drawText(textRect, Qt::AlignRight|Qt::AlignVCenter, QString::number(dB));
        }
    }

    // Draw demod filter box
    if (m_FilterBoxEnabled)
    {
        m_DemodFreqX = xFromFreq(m_DemodCenterFreq);
        m_DemodLowCutFreqX = xFromFreq(m_DemodCenterFreq + m_DemodLowCutFreq);
        m_DemodHiCutFreqX = xFromFreq(m_DemodCenterFreq + m_DemodHiCutFreq);

        int dw = m_DemodHiCutFreqX - m_DemodLowCutFreqX;

        painter.fillRect(m_DemodLowCutFreqX, 0, dw, h,
                         QColor::fromRgba(PLOTTER_FILTER_BOX_COLOR));

        painter.setPen(QPen(QColor::fromRgba(PLOTTER_FILTER_LINE_COLOR), m_DPR));
        painter.drawLine(m_DemodFreqX, 0, m_DemodFreqX, h);
    }

    // Draw a black line at the bottom of the plotter to separate it from the
    // waterfall
    painter.fillRect(QRectF(0.0, h - 1.0 * m_DPR, w, 1.0 * m_DPR), Qt::black);

    painter.end();
}

// Create frequency division strings based on start frequency, span frequency,
// and frequency units.
// Places in QString array m_HDivText
// Keeps all strings the same fractional length
void CPlotter::makeFrequencyStrs()
{
    qint64  StartFreq = m_StartFreqAdj;
    double  freq;
    int     i,j;

    if ((1 == m_FreqUnits) || (m_FreqDigits == 0))
    {
        // if units is Hz then just output integer freq
        for (i = 0; i <= m_HorDivs; i++)
        {
            freq = (double)StartFreq/(double)m_FreqUnits;
            m_HDivText[i].setNum((int)freq);
            StartFreq += m_FreqPerDiv;
        }
        return;
    }
    // here if is fractional frequency values
    // so create max sized text based on frequency units
    for (i = 0; i <= m_HorDivs; i++)
    {
        freq = (double)StartFreq / (double)m_FreqUnits;
        m_HDivText[i].setNum(freq,'f', m_FreqDigits);
        StartFreq += m_FreqPerDiv;
    }
    // now find the division text with the longest non-zero digit
    // to the right of the decimal point.
    int max = 0;
    for (i = 0; i <= m_HorDivs; i++)
    {
        int dp = m_HDivText[i].indexOf('.');
        int l = m_HDivText[i].length()-1;
        for (j = l; j > dp; j--)
        {
            if (m_HDivText[i][j] != '0')
                break;
        }
        if ((j - dp) > max)
            max = j - dp;
    }
    // truncate all strings to maximum fractional length
    StartFreq = m_StartFreqAdj;
    for (i = 0; i <= m_HorDivs; i++)
    {
        freq = (double)StartFreq/(double)m_FreqUnits;
        m_HDivText[i].setNum(freq,'f', max);
        StartFreq += m_FreqPerDiv;
    }
}

// Convert from frequency to screen coordinate
int CPlotter::xFromFreq(qint64 freq)
{
    qreal w = m_Size.width() * m_DPR;
    double startFreq = (double)m_CenterFreq
                       + (double)m_FftCenter
                       - (double)m_Span / 2.0;
    int x = qRound(w * ((double)freq - startFreq) / (double)m_Span);
    return x;
}

// Convert from screen coordinate to frequency
qint64 CPlotter::freqFromX(int x)
{
    double ratio = 0;
    if ((m_Size.width() > 0) && (m_DPR > 0))
        ratio = (double)x / (qreal)m_Size.width() / m_DPR;

    qint64 f = qRound64((double)m_CenterFreq + (double)m_FftCenter
                        - (double)m_Span / 2.0 + ratio * (double)m_Span);
    return f;
}

/** Calculate time offset of a given line on the waterfall */
quint64 CPlotter::msecFromY(int y)
{
    int h = m_OverlayPixmap.height();

    // ensure we are in the waterfall region
    if (y < h)
        return 0;

    qreal dy = (qreal)y - (qreal)h;

    if (msec_per_wfline > 0)
        return tlast_wf_drawn_ms - dy * msec_per_wfline;
    else
        return tlast_wf_drawn_ms - dy * getWfTimeRes();
}

// Round frequency to click resolution value
qint64 CPlotter::roundFreq(qint64 freq, int resolution)
{
    qint64 delta = resolution;
    qint64 delta_2 = delta / 2;
    if (freq >= 0)
        return (freq - (freq + delta_2) % delta + delta_2);
    else
        return (freq - (freq + delta_2) % delta - delta_2);
}

// Clamp demod freqeuency limits of m_DemodCenterFreq
void CPlotter::clampDemodParameters()
{
    if(m_DemodLowCutFreq < m_FLowCmin)
        m_DemodLowCutFreq = m_FLowCmin;
    if(m_DemodLowCutFreq > m_FLowCmax)
        m_DemodLowCutFreq = m_FLowCmax;

    if(m_DemodHiCutFreq < m_FHiCmin)
        m_DemodHiCutFreq = m_FHiCmin;
    if(m_DemodHiCutFreq > m_FHiCmax)
        m_DemodHiCutFreq = m_FHiCmax;
}

void CPlotter::setDemodRanges(int FLowCmin, int FLowCmax,
                              int FHiCmin, int FHiCmax,
                              bool symetric)
{
    m_FLowCmin=FLowCmin;
    m_FLowCmax=FLowCmax;
    m_FHiCmin=FHiCmin;
    m_FHiCmax=FHiCmax;
    m_symetric=symetric;
    clampDemodParameters();
    updateOverlay();
}

void CPlotter::setCenterFreq(quint64 f)
{
    if((quint64)m_CenterFreq == f)
        return;

    qint64 offset = m_CenterFreq - m_DemodCenterFreq;

    m_CenterFreq = f;
    m_DemodCenterFreq = m_CenterFreq - offset;

    m_MaxHoldValid = false;
    m_MinHoldValid = false;
    m_histIIRValid = false;
    m_IIRValid = false;

    updateOverlay();
}

// Invalidate overlay. If not running, force a redraw.
void CPlotter::updateOverlay()
{
    m_DrawOverlay = true;
    draw(false);
}

/** Reset horizontal zoom to 100% and centered around 0. */
void CPlotter::resetHorizontalZoom(void)
{
    setFftCenterFreq(0);
    setSpanFreq((qint32)m_SampleFreq);
    emit newZoomLevel(1.0);
    m_MaxHoldValid = false;
    m_MinHoldValid = false;
    m_histIIRValid = false;
    updateOverlay();
}

/** Center FFT plot around 0 (corresponds to center freq). */
void CPlotter::moveToCenterFreq()
{
    setFftCenterFreq(0);
    m_MaxHoldValid = false;
    m_MinHoldValid = false;
    m_histIIRValid = false;
    updateOverlay();
}

/** Center FFT plot around the demodulator frequency. */
void CPlotter::moveToDemodFreq()
{
    setFftCenterFreq(m_DemodCenterFreq-m_CenterFreq);
    m_MaxHoldValid = false;
    m_MinHoldValid = false;
    m_histIIRValid = false;
    updateOverlay();
}

/** Set FFT plot color. */
void CPlotter::setFftPlotColor(const QColor& color)
{
    m_avgFftColor = color;
    m_maxFftColor = color;
    // m_maxFftColor.setAlpha(192);
    m_FftFillCol = color;
    m_FftFillCol.setAlpha(26);
    m_MaxHoldColor = color;
    m_MaxHoldColor.setAlpha(80);
    m_MinHoldColor = color;
    m_MinHoldColor.setAlpha(80);
}

/** Enable/disable filling the area below the FFT plot. */
void CPlotter::enableFftFill(bool enabled)
{
    m_FftFill = enabled;
}

/** Set peak hold on or off. */
void CPlotter::enableMaxHold(bool enabled)
{
    m_MaxHoldActive = enabled;
    m_MaxHoldValid = false;
}

/** Set min hold on or off. */
void CPlotter::enableMinHold(bool enabled)
{
    m_MinHoldActive = enabled;
    m_MinHoldValid = false;
}

/**
 * Set peak detection on or off.
 * @param enabled The new state of peak detection.
 * @param c Minimum distance of peaks from mean, in multiples of standard deviation.
 */
void CPlotter::enablePeakDetect(bool enabled)
{
    m_PeakDetectActive = enabled;
}

void CPlotter::enableBandPlan(bool enabled)
{
    m_BandPlanEnabled = enabled;
    updateOverlay();
}

void CPlotter::enableMarkers(bool enabled)
{
    m_MarkersEnabled = enabled;
}

void CPlotter::setMarkers(qint64 a, qint64 b)
{
    // Invalidate x positions
    m_MarkerAX = -1;
    m_MarkerBX = -1;

    m_MarkerFreqA = a;
    m_MarkerFreqB = b;

    updateOverlay();
}

void CPlotter::clearWaterfall()
{
    if (!m_WaterfallImage.isNull()) {
        m_WaterfallImage.fill(Qt::black);
    }
}

void CPlotter::calcDivSize (qint64 low, qint64 high, int divswanted, qint64 &adjlow, qint64 &step, int& divs)
{
    qCDebug(plotter) << "low:" << low;
    qCDebug(plotter) << "high:" << high;
    qCDebug(plotter) << "divswanted:" << divswanted;

    if (divswanted == 0)
        return;

    static const qint64 stepTable[] = { 1, 2, 5 };
    static const int stepTableSize = sizeof (stepTable) / sizeof (stepTable[0]);
    qint64 multiplier = 1;
    step = 1;
    divs = high - low;
    int index = 0;
    adjlow = (low / step) * step;

    while (divs > divswanted)
    {
        step = stepTable[index] * multiplier;
        divs = int ((high - low) / step);
        adjlow = (low / step) * step;
        index = index + 1;
        if (index == stepTableSize)
        {
            index = 0;
            multiplier = multiplier * 10;
        }
    }
    if (adjlow < low)
        adjlow += step;

    qCDebug(plotter) << "adjlow:" << adjlow;
    qCDebug(plotter) << "step:" << step;
    qCDebug(plotter) << "divs:" << divs;
}

void CPlotter::showToolTip(QMouseEvent* event, QString toolTipText)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QToolTip::showText(event->globalPos(), toolTipText, this);
#else
    QToolTip::showText(event->globalPosition().toPoint(), toolTipText, this);
#endif
}

// contributed by Chris Kuethe @ckuethe
// source https://ai.googleblog.com/2019/08/turbo-improved-rainbow-colormap-for.html
unsigned char turbo[256][3] = {
        {48,18,59}, {50,21,67},    {51,24,74},
        {52,27,81},   {53,30,88},   {54,33,95},    {55,36,102},   {56,39,109},
        {57,42,115},  {58,45,121},  {59,47,128},   {60,50,134},   {61,53,139},
        {62,56,145},  {63,59,151},  {63,62,156},   {64,64,162},   {65,67,167},
        {65,70,172},  {66,73,177},  {66,75,181},   {67,78,186},   {68,81,191},
        {68,84,195},  {68,86,199},  {69,89,203},   {69,92,207},   {69,94,211},
        {70,97,214},  {70,100,218}, {70,102,221},  {70,105,224},  {70,107,227},
        {71,110,230}, {71,113,233}, {71,115,235},  {71,118,238},  {71,120,240},
        {71,123,242}, {70,125,244}, {70,128,246},  {70,130,248},  {70,133,250},
        {70,135,251}, {69,138,252}, {69,140,253},  {68,143,254},  {67,145,254},
        {66,148,255}, {65,150,255}, {64,153,255},  {62,155,254},  {61,158,254},
        {59,160,253}, {58,163,252}, {56,165,251},  {55,168,250},  {53,171,248},
        {51,173,247}, {49,175,245}, {47,178,244},  {46,180,242},  {44,183,240},
        {42,185,238}, {40,188,235}, {39,190,233},  {37,192,231},  {35,195,228},
        {34,197,226}, {32,199,223}, {31,201,221},  {30,203,218},  {28,205,216},
        {27,208,213}, {26,210,210}, {26,212,208},  {25,213,205},  {24,215,202},
        {24,217,200}, {24,219,197}, {24,221,194},  {24,222,192},  {24,224,189},
        {25,226,187}, {25,227,185}, {26,228,182},  {28,230,180},  {29,231,178},
        {31,233,175}, {32,234,172}, {34,235,170},  {37,236,167},  {39,238,164},
        {42,239,161}, {44,240,158}, {47,241,155},  {50,242,152},  {53,243,148},
        {56,244,145}, {60,245,142}, {63,246,138},  {67,247,135},  {70,248,132},
        {74,248,128}, {78,249,125}, {82,250,122},  {85,250,118},  {89,251,115},
        {93,252,111}, {97,252,108}, {101,253,105}, {105,253,102}, {109,254,98},
        {113,254,95}, {117,254,92}, {121,254,89},  {125,255,86},  {128,255,83},
        {132,255,81}, {136,255,78}, {139,255,75},  {143,255,73},  {146,255,71},
        {150,254,68}, {153,254,66}, {156,254,64},  {159,253,63},  {161,253,61},
        {164,252,60}, {167,252,58}, {169,251,57},  {172,251,56},  {175,250,55},
        {177,249,54}, {180,248,54}, {183,247,53},  {185,246,53},  {188,245,52},
        {190,244,52}, {193,243,52}, {195,241,52},  {198,240,52},  {200,239,52},
        {203,237,52}, {205,236,52}, {208,234,52},  {210,233,53},  {212,231,53},
        {215,229,53}, {217,228,54}, {219,226,54},  {221,224,55},  {223,223,55},
        {225,221,55}, {227,219,56}, {229,217,56},  {231,215,57},  {233,213,57},
        {235,211,57}, {236,209,58}, {238,207,58},  {239,205,58},  {241,203,58},
        {242,201,58}, {244,199,58}, {245,197,58},  {246,195,58},  {247,193,58},
        {248,190,57}, {249,188,57}, {250,186,57},  {251,184,56},  {251,182,55},
        {252,179,54}, {252,177,54}, {253,174,53},  {253,172,52},  {254,169,51},
        {254,167,50}, {254,164,49}, {254,161,48},  {254,158,47},  {254,155,45},
        {254,153,44}, {254,150,43}, {254,147,42},  {254,144,41},  {253,141,39},
        {253,138,38}, {252,135,37}, {252,132,35},  {251,129,34},  {251,126,33},
        {250,123,31}, {249,120,30}, {249,117,29},  {248,114,28},  {247,111,26},
        {246,108,25}, {245,105,24}, {244,102,23},  {243,99,21},   {242,96,20},
        {241,93,19},  {240,91,18},  {239,88,17},   {237,85,16},   {236,83,15},
        {235,80,14},  {234,78,13},  {232,75,12},   {231,73,12},   {229,71,11},
        {228,69,10},  {226,67,10},  {225,65,9},    {223,63,8},    {221,61,8},
        {220,59,7},   {218,57,7},   {216,55,6},    {214,53,6},    {212,51,5},
        {210,49,5},   {208,47,5},   {206,45,4},    {204,43,4},    {202,42,4},
        {200,40,3},   {197,38,3},   {195,37,3},    {193,35,2},    {190,33,2},
        {188,32,2},   {185,30,2},   {183,29,2},    {180,27,1},    {178,26,1},
        {175,24,1},   {172,23,1},   {169,22,1},    {167,20,1},    {164,19,1},
        {161,18,1},   {158,16,1},   {155,15,1},    {152,14,1},    {149,13,1},
        {146,11,1},   {142,10,1},   {139,9,2},     {136,8,2},     {133,7,2},
        {129,6,2},    {126,5,2},    {122,4,3}
};

// contributed by @devnulling
unsigned char plasma[256][3] = {
        {12, 7, 134},   {16, 7, 135},   {19, 6, 137},   {21, 6, 138},   {24, 6, 139},
        {27, 6, 140},   {29, 6, 141},   {31, 5, 142},   {33, 5, 143},   {35, 5, 144},
        {37, 5, 145},   {39, 5, 146},   {41, 5, 147},   {43, 5, 148},   {45, 4, 148},
        {47, 4, 149},   {49, 4, 150},   {51, 4, 151},   {52, 4, 152},   {54, 4, 152},
        {56, 4, 153},   {58, 4, 154},   {59, 3, 154},   {61, 3, 155},   {63, 3, 156},
        {64, 3, 156},   {66, 3, 157},   {68, 3, 158},   {69, 3, 158},   {71, 2, 159},
        {73, 2, 159},   {74, 2, 160},   {76, 2, 161},   {78, 2, 161},   {79, 2, 162},
        {81, 1, 162},   {82, 1, 163},   {84, 1, 163},   {86, 1, 163},   {87, 1, 164},
        {89, 1, 164},   {90, 0, 165},   {92, 0, 165},   {94, 0, 165},   {95, 0, 166},
        {97, 0, 166},   {98, 0, 166},   {100, 0, 167},  {101, 0, 167},  {103, 0, 167},
        {104, 0, 167},  {106, 0, 167},  {108, 0, 168},  {109, 0, 168},  {111, 0, 168},
        {112, 0, 168},  {114, 0, 168},  {115, 0, 168},  {117, 0, 168},  {118, 1, 168},
        {120, 1, 168},  {121, 1, 168},  {123, 2, 168},  {124, 2, 167},  {126, 3, 167},
        {127, 3, 167},  {129, 4, 167},  {130, 4, 167},  {132, 5, 166},  {133, 6, 166},
        {134, 7, 166},  {136, 7, 165},  {137, 8, 165},  {139, 9, 164},  {140, 10, 164},
        {142, 12, 164}, {143, 13, 163}, {144, 14, 163}, {146, 15, 162}, {147, 16, 161},
        {149, 17, 161}, {150, 18, 160}, {151, 19, 160}, {153, 20, 159}, {154, 21, 158},
        {155, 23, 158}, {157, 24, 157}, {158, 25, 156}, {159, 26, 155}, {160, 27, 155},
        {162, 28, 154}, {163, 29, 153}, {164, 30, 152}, {165, 31, 151}, {167, 33, 151},
        {168, 34, 150}, {169, 35, 149}, {170, 36, 148}, {172, 37, 147}, {173, 38, 146},
        {174, 39, 145}, {175, 40, 144}, {176, 42, 143}, {177, 43, 143}, {178, 44, 142},
        {180, 45, 141}, {181, 46, 140}, {182, 47, 139}, {183, 48, 138}, {184, 50, 137},
        {185, 51, 136}, {186, 52, 135}, {187, 53, 134}, {188, 54, 133}, {189, 55, 132},
        {190, 56, 131}, {191, 57, 130}, {192, 59, 129}, {193, 60, 128}, {194, 61, 128},
        {195, 62, 127}, {196, 63, 126}, {197, 64, 125}, {198, 65, 124}, {199, 66, 123},
        {200, 68, 122}, {201, 69, 121}, {202, 70, 120}, {203, 71, 119}, {204, 72, 118},
        {205, 73, 117}, {206, 74, 117}, {207, 75, 116}, {208, 77, 115}, {209, 78, 114},
        {209, 79, 113}, {210, 80, 112}, {211, 81, 111}, {212, 82, 110}, {213, 83, 109},
        {214, 85, 109}, {215, 86, 108}, {215, 87, 107}, {216, 88, 106}, {217, 89, 105},
        {218, 90, 104}, {219, 91, 103}, {220, 93, 102}, {220, 94, 102}, {221, 95, 101},
        {222, 96, 100}, {223, 97, 99},  {223, 98, 98},  {224, 100, 97}, {225, 101, 96},
        {226, 102, 96}, {227, 103, 95}, {227, 104, 94}, {228, 106, 93}, {229, 107, 92},
        {229, 108, 91}, {230, 109, 90}, {231, 110, 90}, {232, 112, 89}, {232, 113, 88},
        {233, 114, 87}, {234, 115, 86}, {234, 116, 85}, {235, 118, 84}, {236, 119, 84},
        {236, 120, 83}, {237, 121, 82}, {237, 123, 81}, {238, 124, 80}, {239, 125, 79},
        {239, 126, 78}, {240, 128, 77}, {240, 129, 77}, {241, 130, 76}, {242, 132, 75},
        {242, 133, 74}, {243, 134, 73}, {243, 135, 72}, {244, 137, 71}, {244, 138, 71},
        {245, 139, 70}, {245, 141, 69}, {246, 142, 68}, {246, 143, 67}, {246, 145, 66},
        {247, 146, 65}, {247, 147, 65}, {248, 149, 64}, {248, 150, 63}, {248, 152, 62},
        {249, 153, 61}, {249, 154, 60}, {250, 156, 59}, {250, 157, 58}, {250, 159, 58},
        {250, 160, 57}, {251, 162, 56}, {251, 163, 55}, {251, 164, 54}, {252, 166, 53},
        {252, 167, 53}, {252, 169, 52}, {252, 170, 51}, {252, 172, 50}, {252, 173, 49},
        {253, 175, 49}, {253, 176, 48}, {253, 178, 47}, {253, 179, 46}, {253, 181, 45},
        {253, 182, 45}, {253, 184, 44}, {253, 185, 43}, {253, 187, 43}, {253, 188, 42},
        {253, 190, 41}, {253, 192, 41}, {253, 193, 40}, {253, 195, 40}, {253, 196, 39},
        {253, 198, 38}, {252, 199, 38}, {252, 201, 38}, {252, 203, 37}, {252, 204, 37},
        {252, 206, 37}, {251, 208, 36}, {251, 209, 36}, {251, 211, 36}, {250, 213, 36},
        {250, 214, 36}, {250, 216, 36}, {249, 217, 36}, {249, 219, 36}, {248, 221, 36},
        {248, 223, 36}, {247, 224, 36}, {247, 226, 37}, {246, 228, 37}, {246, 229, 37},
        {245, 231, 38}, {245, 233, 38}, {244, 234, 38}, {243, 236, 38}, {243, 238, 38},
        {242, 240, 38}, {242, 241, 38}, {241, 243, 38}, {240, 245, 37}, {240, 246, 35},
        {239, 248, 33}
};

// contributed by @Piruzzolo
float viridis[256][3] = {
        { 0.267004, 0.004874, 0.329415 },
        { 0.268510, 0.009605, 0.335427 },
        { 0.269944, 0.014625, 0.341379 },
        { 0.271305, 0.019942, 0.347269 },
        { 0.272594, 0.025563, 0.353093 },
        { 0.273809, 0.031497, 0.358853 },
        { 0.274952, 0.037752, 0.364543 },
        { 0.276022, 0.044167, 0.370164 },
        { 0.277018, 0.050344, 0.375715 },
        { 0.277941, 0.056324, 0.381191 },
        { 0.278791, 0.062145, 0.386592 },
        { 0.279566, 0.067836, 0.391917 },
        { 0.280267, 0.073417, 0.397163 },
        { 0.280894, 0.078907, 0.402329 },
        { 0.281446, 0.084320, 0.407414 },
        { 0.281924, 0.089666, 0.412415 },
        { 0.282327, 0.094955, 0.417331 },
        { 0.282656, 0.100196, 0.422160 },
        { 0.282910, 0.105393, 0.426902 },
        { 0.283091, 0.110553, 0.431554 },
        { 0.283197, 0.115680, 0.436115 },
        { 0.283229, 0.120777, 0.440584 },
        { 0.283187, 0.125848, 0.444960 },
        { 0.283072, 0.130895, 0.449241 },
        { 0.282884, 0.135920, 0.453427 },
        { 0.282623, 0.140926, 0.457517 },
        { 0.282290, 0.145912, 0.461510 },
        { 0.281887, 0.150881, 0.465405 },
        { 0.281412, 0.155834, 0.469201 },
        { 0.280868, 0.160771, 0.472899 },
        { 0.280255, 0.165693, 0.476498 },
        { 0.279574, 0.170599, 0.479997 },
        { 0.278826, 0.175490, 0.483397 },
        { 0.278012, 0.180367, 0.486697 },
        { 0.277134, 0.185228, 0.489898 },
        { 0.276194, 0.190074, 0.493001 },
        { 0.275191, 0.194905, 0.496005 },
        { 0.274128, 0.199721, 0.498911 },
        { 0.273006, 0.204520, 0.501721 },
        { 0.271828, 0.209303, 0.504434 },
        { 0.270595, 0.214069, 0.507052 },
        { 0.269308, 0.218818, 0.509577 },
        { 0.267968, 0.223549, 0.512008 },
        { 0.266580, 0.228262, 0.514349 },
        { 0.265145, 0.232956, 0.516599 },
        { 0.263663, 0.237631, 0.518762 },
        { 0.262138, 0.242286, 0.520837 },
        { 0.260571, 0.246922, 0.522828 },
        { 0.258965, 0.251537, 0.524736 },
        { 0.257322, 0.256130, 0.526563 },
        { 0.255645, 0.260703, 0.528312 },
        { 0.253935, 0.265254, 0.529983 },
        { 0.252194, 0.269783, 0.531579 },
        { 0.250425, 0.274290, 0.533103 },
        { 0.248629, 0.278775, 0.534556 },
        { 0.246811, 0.283237, 0.535941 },
        { 0.244972, 0.287675, 0.537260 },
        { 0.243113, 0.292092, 0.538516 },
        { 0.241237, 0.296485, 0.539709 },
        { 0.239346, 0.300855, 0.540844 },
        { 0.237441, 0.305202, 0.541921 },
        { 0.235526, 0.309527, 0.542944 },
        { 0.233603, 0.313828, 0.543914 },
        { 0.231674, 0.318106, 0.544834 },
        { 0.229739, 0.322361, 0.545706 },
        { 0.227802, 0.326594, 0.546532 },
        { 0.225863, 0.330805, 0.547314 },
        { 0.223925, 0.334994, 0.548053 },
        { 0.221989, 0.339161, 0.548752 },
        { 0.220057, 0.343307, 0.549413 },
        { 0.218130, 0.347432, 0.550038 },
        { 0.216210, 0.351535, 0.550627 },
        { 0.214298, 0.355619, 0.551184 },
        { 0.212395, 0.359683, 0.551710 },
        { 0.210503, 0.363727, 0.552206 },
        { 0.208623, 0.367752, 0.552675 },
        { 0.206756, 0.371758, 0.553117 },
        { 0.204903, 0.375746, 0.553533 },
        { 0.203063, 0.379716, 0.553925 },
        { 0.201239, 0.383670, 0.554294 },
        { 0.199430, 0.387607, 0.554642 },
        { 0.197636, 0.391528, 0.554969 },
        { 0.195860, 0.395433, 0.555276 },
        { 0.194100, 0.399323, 0.555565 },
        { 0.192357, 0.403199, 0.555836 },
        { 0.190631, 0.407061, 0.556089 },
        { 0.188923, 0.410910, 0.556326 },
        { 0.187231, 0.414746, 0.556547 },
        { 0.185556, 0.418570, 0.556753 },
        { 0.183898, 0.422383, 0.556944 },
        { 0.182256, 0.426184, 0.557120 },
        { 0.180629, 0.429975, 0.557282 },
        { 0.179019, 0.433756, 0.557430 },
        { 0.177423, 0.437527, 0.557565 },
        { 0.175841, 0.441290, 0.557685 },
        { 0.174274, 0.445044, 0.557792 },
        { 0.172719, 0.448791, 0.557885 },
        { 0.171176, 0.452530, 0.557965 },
        { 0.169646, 0.456262, 0.558030 },
        { 0.168126, 0.459988, 0.558082 },
        { 0.166617, 0.463708, 0.558119 },
        { 0.165117, 0.467423, 0.558141 },
        { 0.163625, 0.471133, 0.558148 },
        { 0.162142, 0.474838, 0.558140 },
        { 0.160665, 0.478540, 0.558115 },
        { 0.159194, 0.482237, 0.558073 },
        { 0.157729, 0.485932, 0.558013 },
        { 0.156270, 0.489624, 0.557936 },
        { 0.154815, 0.493313, 0.557840 },
        { 0.153364, 0.497000, 0.557724 },
        { 0.151918, 0.500685, 0.557587 },
        { 0.150476, 0.504369, 0.557430 },
        { 0.149039, 0.508051, 0.557250 },
        { 0.147607, 0.511733, 0.557049 },
        { 0.146180, 0.515413, 0.556823 },
        { 0.144759, 0.519093, 0.556572 },
        { 0.143343, 0.522773, 0.556295 },
        { 0.141935, 0.526453, 0.555991 },
        { 0.140536, 0.530132, 0.555659 },
        { 0.139147, 0.533812, 0.555298 },
        { 0.137770, 0.537492, 0.554906 },
        { 0.136408, 0.541173, 0.554483 },
        { 0.135066, 0.544853, 0.554029 },
        { 0.133743, 0.548535, 0.553541 },
        { 0.132444, 0.552216, 0.553018 },
        { 0.131172, 0.555899, 0.552459 },
        { 0.129933, 0.559582, 0.551864 },
        { 0.128729, 0.563265, 0.551229 },
        { 0.127568, 0.566949, 0.550556 },
        { 0.126453, 0.570633, 0.549841 },
        { 0.125394, 0.574318, 0.549086 },
        { 0.124395, 0.578002, 0.548287 },
        { 0.123463, 0.581687, 0.547445 },
        { 0.122606, 0.585371, 0.546557 },
        { 0.121831, 0.589055, 0.545623 },
        { 0.121148, 0.592739, 0.544641 },
        { 0.120565, 0.596422, 0.543611 },
        { 0.120092, 0.600104, 0.542530 },
        { 0.119738, 0.603785, 0.541400 },
        { 0.119512, 0.607464, 0.540218 },
        { 0.119423, 0.611141, 0.538982 },
        { 0.119483, 0.614817, 0.537692 },
        { 0.119699, 0.618490, 0.536347 },
        { 0.120081, 0.622161, 0.534946 },
        { 0.120638, 0.625828, 0.533488 },
        { 0.121380, 0.629492, 0.531973 },
        { 0.122312, 0.633153, 0.530398 },
        { 0.123444, 0.636809, 0.528763 },
        { 0.124780, 0.640461, 0.527068 },
        { 0.126326, 0.644107, 0.525311 },
        { 0.128087, 0.647749, 0.523491 },
        { 0.130067, 0.651384, 0.521608 },
        { 0.132268, 0.655014, 0.519661 },
        { 0.134692, 0.658636, 0.517649 },
        { 0.137339, 0.662252, 0.515571 },
        { 0.140210, 0.665859, 0.513427 },
        { 0.143303, 0.669459, 0.511215 },
        { 0.146616, 0.673050, 0.508936 },
        { 0.150148, 0.676631, 0.506589 },
        { 0.153894, 0.680203, 0.504172 },
        { 0.157851, 0.683765, 0.501686 },
        { 0.162016, 0.687316, 0.499129 },
        { 0.166383, 0.690856, 0.496502 },
        { 0.170948, 0.694384, 0.493803 },
        { 0.175707, 0.697900, 0.491033 },
        { 0.180653, 0.701402, 0.488189 },
        { 0.185783, 0.704891, 0.485273 },
        { 0.191090, 0.708366, 0.482284 },
        { 0.196571, 0.711827, 0.479221 },
        { 0.202219, 0.715272, 0.476084 },
        { 0.208030, 0.718701, 0.472873 },
        { 0.214000, 0.722114, 0.469588 },
        { 0.220124, 0.725509, 0.466226 },
        { 0.226397, 0.728888, 0.462789 },
        { 0.232815, 0.732247, 0.459277 },
        { 0.239374, 0.735588, 0.455688 },
        { 0.246070, 0.738910, 0.452024 },
        { 0.252899, 0.742211, 0.448284 },
        { 0.259857, 0.745492, 0.444467 },
        { 0.266941, 0.748751, 0.440573 },
        { 0.274149, 0.751988, 0.436601 },
        { 0.281477, 0.755203, 0.432552 },
        { 0.288921, 0.758394, 0.428426 },
        { 0.296479, 0.761561, 0.424223 },
        { 0.304148, 0.764704, 0.419943 },
        { 0.311925, 0.767822, 0.415586 },
        { 0.319809, 0.770914, 0.411152 },
        { 0.327796, 0.773980, 0.406640 },
        { 0.335885, 0.777018, 0.402049 },
        { 0.344074, 0.780029, 0.397381 },
        { 0.352360, 0.783011, 0.392636 },
        { 0.360741, 0.785964, 0.387814 },
        { 0.369214, 0.788888, 0.382914 },
        { 0.377779, 0.791781, 0.377939 },
        { 0.386433, 0.794644, 0.372886 },
        { 0.395174, 0.797475, 0.367757 },
        { 0.404001, 0.800275, 0.362552 },
        { 0.412913, 0.803041, 0.357269 },
        { 0.421908, 0.805774, 0.351910 },
        { 0.430983, 0.808473, 0.346476 },
        { 0.440137, 0.811138, 0.340967 },
        { 0.449368, 0.813768, 0.335384 },
        { 0.458674, 0.816363, 0.329727 },
        { 0.468053, 0.818921, 0.323998 },
        { 0.477504, 0.821444, 0.318195 },
        { 0.487026, 0.823929, 0.312321 },
        { 0.496615, 0.826376, 0.306377 },
        { 0.506271, 0.828786, 0.300362 },
        { 0.515992, 0.831158, 0.294279 },
        { 0.525776, 0.833491, 0.288127 },
        { 0.535621, 0.835785, 0.281908 },
        { 0.545524, 0.838039, 0.275626 },
        { 0.555484, 0.840254, 0.269281 },
        { 0.565498, 0.842430, 0.262877 },
        { 0.575563, 0.844566, 0.256415 },
        { 0.585678, 0.846661, 0.249897 },
        { 0.595839, 0.848717, 0.243329 },
        { 0.606045, 0.850733, 0.236712 },
        { 0.616293, 0.852709, 0.230052 },
        { 0.626579, 0.854645, 0.223353 },
        { 0.636902, 0.856542, 0.216620 },
        { 0.647257, 0.858400, 0.209861 },
        { 0.657642, 0.860219, 0.203082 },
        { 0.668054, 0.861999, 0.196293 },
        { 0.678489, 0.863742, 0.189503 },
        { 0.688944, 0.865448, 0.182725 },
        { 0.699415, 0.867117, 0.175971 },
        { 0.709898, 0.868751, 0.169257 },
        { 0.720391, 0.870350, 0.162603 },
        { 0.730889, 0.871916, 0.156029 },
        { 0.741388, 0.873449, 0.149561 },
        { 0.751884, 0.874951, 0.143228 },
        { 0.762373, 0.876424, 0.137064 },
        { 0.772852, 0.877868, 0.131109 },
        { 0.783315, 0.879285, 0.125405 },
        { 0.793760, 0.880678, 0.120005 },
        { 0.804182, 0.882046, 0.114965 },
        { 0.814576, 0.883393, 0.110347 },
        { 0.824940, 0.884720, 0.106217 },
        { 0.835270, 0.886029, 0.102646 },
        { 0.845561, 0.887322, 0.099702 },
        { 0.855810, 0.888601, 0.097452 },
        { 0.866013, 0.889868, 0.095953 },
        { 0.876168, 0.891125, 0.095250 },
        { 0.886271, 0.892374, 0.095374 },
        { 0.896320, 0.893616, 0.096335 },
        { 0.906311, 0.894855, 0.098125 },
        { 0.916242, 0.896091, 0.100717 },
        { 0.926106, 0.897330, 0.104071 },
        { 0.935904, 0.898570, 0.108131 },
        { 0.945636, 0.899815, 0.112838 },
        { 0.955300, 0.901065, 0.118128 },
        { 0.964894, 0.902323, 0.123941 },
        { 0.974417, 0.903590, 0.130215 },
        { 0.983868, 0.904867, 0.136897 },
        { 0.993248, 0.906157, 0.143936 }
};

// contributed by @cmj
float magma[256][3] = {
	{ 0.001462, 0.000466, 0.013866 },
	{ 0.002258, 0.001295, 0.018331 },
	{ 0.003279, 0.002305, 0.023708 },
	{ 0.004512, 0.003490, 0.029965 },
	{ 0.005950, 0.004843, 0.037130 },
	{ 0.007588, 0.006356, 0.044973 },
	{ 0.009426, 0.008022, 0.052844 },
	{ 0.011465, 0.009828, 0.060750 },
	{ 0.013708, 0.011771, 0.068667 },
	{ 0.016156, 0.013840, 0.076603 },
	{ 0.018815, 0.016026, 0.084584 },
	{ 0.021692, 0.018320, 0.092610 },
	{ 0.024792, 0.020715, 0.100676 },
	{ 0.028123, 0.023201, 0.108787 },
	{ 0.031696, 0.025765, 0.116965 },
	{ 0.035520, 0.028397, 0.125209 },
	{ 0.039608, 0.031090, 0.133515 },
	{ 0.043830, 0.033830, 0.141886 },
	{ 0.048062, 0.036607, 0.150327 },
	{ 0.052320, 0.039407, 0.158841 },
	{ 0.056615, 0.042160, 0.167446 },
	{ 0.060949, 0.044794, 0.176129 },
	{ 0.065330, 0.047318, 0.184892 },
	{ 0.069764, 0.049726, 0.193735 },
	{ 0.074257, 0.052017, 0.202660 },
	{ 0.078815, 0.054184, 0.211667 },
	{ 0.083446, 0.056225, 0.220755 },
	{ 0.088155, 0.058133, 0.229922 },
	{ 0.092949, 0.059904, 0.239164 },
	{ 0.097833, 0.061531, 0.248477 },
	{ 0.102815, 0.063010, 0.257854 },
	{ 0.107899, 0.064335, 0.267289 },
	{ 0.113094, 0.065492, 0.276784 },
	{ 0.118405, 0.066479, 0.286321 },
	{ 0.123833, 0.067295, 0.295879 },
	{ 0.129380, 0.067935, 0.305443 },
	{ 0.135053, 0.068391, 0.315000 },
	{ 0.140858, 0.068654, 0.324538 },
	{ 0.146785, 0.068738, 0.334011 },
	{ 0.152839, 0.068637, 0.343404 },
	{ 0.159018, 0.068354, 0.352688 },
	{ 0.165308, 0.067911, 0.361816 },
	{ 0.171713, 0.067305, 0.370771 },
	{ 0.178212, 0.066576, 0.379497 },
	{ 0.184801, 0.065732, 0.387973 },
	{ 0.191460, 0.064818, 0.396152 },
	{ 0.198177, 0.063862, 0.404009 },
	{ 0.204935, 0.062907, 0.411514 },
	{ 0.211718, 0.061992, 0.418647 },
	{ 0.218512, 0.061158, 0.425392 },
	{ 0.225302, 0.060445, 0.431742 },
	{ 0.232077, 0.059889, 0.437695 },
	{ 0.238826, 0.059517, 0.443256 },
	{ 0.245543, 0.059352, 0.448436 },
	{ 0.252220, 0.059415, 0.453248 },
	{ 0.258857, 0.059706, 0.457710 },
	{ 0.265447, 0.060237, 0.461840 },
	{ 0.271994, 0.060994, 0.465660 },
	{ 0.278493, 0.061978, 0.469190 },
	{ 0.284951, 0.063168, 0.472451 },
	{ 0.291366, 0.064553, 0.475462 },
	{ 0.297740, 0.066117, 0.478243 },
	{ 0.304081, 0.067835, 0.480812 },
	{ 0.310382, 0.069702, 0.483186 },
	{ 0.316654, 0.071690, 0.485380 },
	{ 0.322899, 0.073782, 0.487408 },
	{ 0.329114, 0.075972, 0.489287 },
	{ 0.335308, 0.078236, 0.491024 },
	{ 0.341482, 0.080564, 0.492631 },
	{ 0.347636, 0.082946, 0.494121 },
	{ 0.353773, 0.085373, 0.495501 },
	{ 0.359898, 0.087831, 0.496778 },
	{ 0.366012, 0.090314, 0.497960 },
	{ 0.372116, 0.092816, 0.499053 },
	{ 0.378211, 0.095332, 0.500067 },
	{ 0.384299, 0.097855, 0.501002 },
	{ 0.390384, 0.100379, 0.501864 },
	{ 0.396467, 0.102902, 0.502658 },
	{ 0.402548, 0.105420, 0.503386 },
	{ 0.408629, 0.107930, 0.504052 },
	{ 0.414709, 0.110431, 0.504662 },
	{ 0.420791, 0.112920, 0.505215 },
	{ 0.426877, 0.115395, 0.505714 },
	{ 0.432967, 0.117855, 0.506160 },
	{ 0.439062, 0.120298, 0.506555 },
	{ 0.445163, 0.122724, 0.506901 },
	{ 0.451271, 0.125132, 0.507198 },
	{ 0.457386, 0.127522, 0.507448 },
	{ 0.463508, 0.129893, 0.507652 },
	{ 0.469640, 0.132245, 0.507809 },
	{ 0.475780, 0.134577, 0.507921 },
	{ 0.481929, 0.136891, 0.507989 },
	{ 0.488088, 0.139186, 0.508011 },
	{ 0.494258, 0.141462, 0.507988 },
	{ 0.500438, 0.143719, 0.507920 },
	{ 0.506629, 0.145958, 0.507806 },
	{ 0.512831, 0.148179, 0.507648 },
	{ 0.519045, 0.150383, 0.507443 },
	{ 0.525270, 0.152569, 0.507192 },
	{ 0.531507, 0.154739, 0.506895 },
	{ 0.537755, 0.156894, 0.506551 },
	{ 0.544015, 0.159033, 0.506159 },
	{ 0.550287, 0.161158, 0.505719 },
	{ 0.556571, 0.163269, 0.505230 },
	{ 0.562866, 0.165368, 0.504692 },
	{ 0.569172, 0.167454, 0.504105 },
	{ 0.575490, 0.169530, 0.503466 },
	{ 0.581819, 0.171596, 0.502777 },
	{ 0.588158, 0.173652, 0.502035 },
	{ 0.594508, 0.175701, 0.501241 },
	{ 0.600868, 0.177743, 0.500394 },
	{ 0.607238, 0.179779, 0.499492 },
	{ 0.613617, 0.181811, 0.498536 },
	{ 0.620005, 0.183840, 0.497524 },
	{ 0.626401, 0.185867, 0.496456 },
	{ 0.632805, 0.187893, 0.495332 },
	{ 0.639216, 0.189921, 0.494150 },
	{ 0.645633, 0.191952, 0.492910 },
	{ 0.652056, 0.193986, 0.491611 },
	{ 0.658483, 0.196027, 0.490253 },
	{ 0.664915, 0.198075, 0.488836 },
	{ 0.671349, 0.200133, 0.487358 },
	{ 0.677786, 0.202203, 0.485819 },
	{ 0.684224, 0.204286, 0.484219 },
	{ 0.690661, 0.206384, 0.482558 },
	{ 0.697098, 0.208501, 0.480835 },
	{ 0.703532, 0.210638, 0.479049 },
	{ 0.709962, 0.212797, 0.477201 },
	{ 0.716387, 0.214982, 0.475290 },
	{ 0.722805, 0.217194, 0.473316 },
	{ 0.729216, 0.219437, 0.471279 },
	{ 0.735616, 0.221713, 0.469180 },
	{ 0.742004, 0.224025, 0.467018 },
	{ 0.748378, 0.226377, 0.464794 },
	{ 0.754737, 0.228772, 0.462509 },
	{ 0.761077, 0.231214, 0.460162 },
	{ 0.767398, 0.233705, 0.457755 },
	{ 0.773695, 0.236249, 0.455289 },
	{ 0.779968, 0.238851, 0.452765 },
	{ 0.786212, 0.241514, 0.450184 },
	{ 0.792427, 0.244242, 0.447543 },
	{ 0.798608, 0.247040, 0.444848 },
	{ 0.804752, 0.249911, 0.442102 },
	{ 0.810855, 0.252861, 0.439305 },
	{ 0.816914, 0.255895, 0.436461 },
	{ 0.822926, 0.259016, 0.433573 },
	{ 0.828886, 0.262229, 0.430644 },
	{ 0.834791, 0.265540, 0.427671 },
	{ 0.840636, 0.268953, 0.424666 },
	{ 0.846416, 0.272473, 0.421631 },
	{ 0.852126, 0.276106, 0.418573 },
	{ 0.857763, 0.279857, 0.415496 },
	{ 0.863320, 0.283729, 0.412403 },
	{ 0.868793, 0.287728, 0.409303 },
	{ 0.874176, 0.291859, 0.406205 },
	{ 0.879464, 0.296125, 0.403118 },
	{ 0.884651, 0.300530, 0.400047 },
	{ 0.889731, 0.305079, 0.397002 },
	{ 0.894700, 0.309773, 0.393995 },
	{ 0.899552, 0.314616, 0.391037 },
	{ 0.904281, 0.319610, 0.388137 },
	{ 0.908884, 0.324755, 0.385308 },
	{ 0.913354, 0.330052, 0.382563 },
	{ 0.917689, 0.335500, 0.379915 },
	{ 0.921884, 0.341098, 0.377376 },
	{ 0.925937, 0.346844, 0.374959 },
	{ 0.929845, 0.352734, 0.372677 },
	{ 0.933606, 0.358764, 0.370541 },
	{ 0.937221, 0.364929, 0.368567 },
	{ 0.940687, 0.371224, 0.366762 },
	{ 0.944006, 0.377643, 0.365136 },
	{ 0.947180, 0.384178, 0.363701 },
	{ 0.950210, 0.390820, 0.362468 },
	{ 0.953099, 0.397563, 0.361438 },
	{ 0.955849, 0.404400, 0.360619 },
	{ 0.958464, 0.411324, 0.360014 },
	{ 0.960949, 0.418323, 0.359630 },
	{ 0.963310, 0.425390, 0.359469 },
	{ 0.965549, 0.432519, 0.359529 },
	{ 0.967671, 0.439703, 0.359810 },
	{ 0.969680, 0.446936, 0.360311 },
	{ 0.971582, 0.454210, 0.361030 },
	{ 0.973381, 0.461520, 0.361965 },
	{ 0.975082, 0.468861, 0.363111 },
	{ 0.976690, 0.476226, 0.364466 },
	{ 0.978210, 0.483612, 0.366025 },
	{ 0.979645, 0.491014, 0.367783 },
	{ 0.981000, 0.498428, 0.369734 },
	{ 0.982279, 0.505851, 0.371874 },
	{ 0.983485, 0.513280, 0.374198 },
	{ 0.984622, 0.520713, 0.376698 },
	{ 0.985693, 0.528148, 0.379371 },
	{ 0.986700, 0.535582, 0.382210 },
	{ 0.987646, 0.543015, 0.385210 },
	{ 0.988533, 0.550446, 0.388365 },
	{ 0.989363, 0.557873, 0.391671 },
	{ 0.990138, 0.565296, 0.395122 },
	{ 0.990871, 0.572706, 0.398714 },
	{ 0.991558, 0.580107, 0.402441 },
	{ 0.992196, 0.587502, 0.406299 },
	{ 0.992785, 0.594891, 0.410283 },
	{ 0.993326, 0.602275, 0.414390 },
	{ 0.993834, 0.609644, 0.418613 },
	{ 0.994309, 0.616999, 0.422950 },
	{ 0.994738, 0.624350, 0.427397 },
	{ 0.995122, 0.631696, 0.431951 },
	{ 0.995480, 0.639027, 0.436607 },
	{ 0.995810, 0.646344, 0.441361 },
	{ 0.996096, 0.653659, 0.446213 },
	{ 0.996341, 0.660969, 0.451160 },
	{ 0.996580, 0.668256, 0.456192 },
	{ 0.996775, 0.675541, 0.461314 },
	{ 0.996925, 0.682828, 0.466526 },
	{ 0.997077, 0.690088, 0.471811 },
	{ 0.997186, 0.697349, 0.477182 },
	{ 0.997254, 0.704611, 0.482635 },
	{ 0.997325, 0.711848, 0.488154 },
	{ 0.997351, 0.719089, 0.493755 },
	{ 0.997351, 0.726324, 0.499428 },
	{ 0.997341, 0.733545, 0.505167 },
	{ 0.997285, 0.740772, 0.510983 },
	{ 0.997228, 0.747981, 0.516859 },
	{ 0.997138, 0.755190, 0.522806 },
	{ 0.997019, 0.762398, 0.528821 },
	{ 0.996898, 0.769591, 0.534892 },
	{ 0.996727, 0.776795, 0.541039 },
	{ 0.996571, 0.783977, 0.547233 },
	{ 0.996369, 0.791167, 0.553499 },
	{ 0.996162, 0.798348, 0.559820 },
	{ 0.995932, 0.805527, 0.566202 },
	{ 0.995680, 0.812706, 0.572645 },
	{ 0.995424, 0.819875, 0.579140 },
	{ 0.995131, 0.827052, 0.585701 },
	{ 0.994851, 0.834213, 0.592307 },
	{ 0.994524, 0.841387, 0.598983 },
	{ 0.994222, 0.848540, 0.605696 },
	{ 0.993866, 0.855711, 0.612482 },
	{ 0.993545, 0.862859, 0.619299 },
	{ 0.993170, 0.870024, 0.626189 },
	{ 0.992831, 0.877168, 0.633109 },
	{ 0.992440, 0.884330, 0.640099 },
	{ 0.992089, 0.891470, 0.647116 },
	{ 0.991688, 0.898627, 0.654202 },
	{ 0.991332, 0.905763, 0.661309 },
	{ 0.990930, 0.912915, 0.668481 },
	{ 0.990570, 0.920049, 0.675675 },
	{ 0.990175, 0.927196, 0.682926 },
	{ 0.989815, 0.934329, 0.690198 },
	{ 0.989434, 0.941470, 0.697519 },
	{ 0.989077, 0.948604, 0.704863 },
	{ 0.988717, 0.955742, 0.712242 },
	{ 0.988367, 0.962878, 0.719649 },
	{ 0.988033, 0.970012, 0.727077 },
	{ 0.987691, 0.977154, 0.734536 },
	{ 0.987387, 0.984288, 0.742002 },
	{ 0.987053, 0.991438, 0.749504 }
};

// contributed by @cmj
float grape[256][3] = {
	{ 0.17842553, 0.13513835, 0.16166147 }, 
	{ 0.18388264, 0.13648638, 0.1646808  }, 
	{ 0.18936728, 0.13778491, 0.16770558 }, 
	{ 0.19488022, 0.13903261, 0.17073898 }, 
	{ 0.2004213 , 0.14022861, 0.17378394 }, 
	{ 0.20598563, 0.14137481, 0.17684278 }, 
	{ 0.21157859, 0.14246709, 0.17991876 }, 
	{ 0.21719224, 0.14350927, 0.18301362 }, 
	{ 0.22283364, 0.144496  , 0.18613078 }, 
	{ 0.22849466, 0.14543139, 0.18927175 }, 
	{ 0.2341828 , 0.1463095 , 0.1924401  }, 
	{ 0.23989079, 0.14713399, 0.19563734 }, 
	{ 0.24562333, 0.14790049, 0.19886677 }, 
	{ 0.25137776, 0.14860955, 0.20213061 }, 
	{ 0.25715275, 0.14926082, 0.20543125 }, 
	{ 0.26295307, 0.14984937, 0.20877243 }, 
	{ 0.26877366, 0.15037739, 0.21215597 }, 
	{ 0.27461474, 0.15084315, 0.21558481 }, 
	{ 0.2804794 , 0.15124248, 0.21906277 }, 
	{ 0.28636528, 0.15157537, 0.22259249 }, 
	{ 0.2922715 , 0.15184054, 0.22617704 }, 
	{ 0.29819806, 0.15203591, 0.22981989 }, 
	{ 0.30414483, 0.15215932, 0.23352464 }, 
	{ 0.31011161, 0.15220851, 0.23729507 }, 
	{ 0.3160985 , 0.15218073, 0.24113528 }, 
	{ 0.32210445, 0.15207411, 0.24504918 }, 
	{ 0.32812853, 0.15188652, 0.24904092 }, 
	{ 0.33417037, 0.15161504, 0.25311512 }, 
	{ 0.34022963, 0.15125649, 0.25727673 }, 
	{ 0.34630479, 0.1508086 , 0.26153044 }, 
	{ 0.35239459, 0.15026865, 0.26588132 }, 
	{ 0.35849751, 0.14963392, 0.27033462 }, 
	{ 0.36461184, 0.14890165, 0.27489576 }, 
	{ 0.3707392 , 0.14806473, 0.27957272 }, 
	{ 0.37687445, 0.14712369, 0.28436966 }, 
	{ 0.38301619, 0.14607434, 0.28929339 }, 
	{ 0.38916394, 0.14491071, 0.29435205 }, 
	{ 0.39531286, 0.14363215, 0.29955093 }, 
	{ 0.40146175, 0.14223269, 0.30489865 }, 
	{ 0.40760722, 0.14070879, 0.31040263 }, 
	{ 0.41374455, 0.13905853, 0.3160695  }, 
	{ 0.41987027, 0.13727763, 0.32190751 }, 
	{ 0.42598098, 0.13536106, 0.32792574 }, 
	{ 0.43207074, 0.13330738, 0.33413104 }, 
	{ 0.43813387, 0.13111453, 0.34053083 }, 
	{ 0.44416424, 0.12878087, 0.3471324  }, 
	{ 0.45015613, 0.12630362, 0.35394403 }, 
	{ 0.45610114, 0.12368516, 0.36097045 }, 
	{ 0.46199101, 0.12092812, 0.36821625 }, 
	{ 0.4678179 , 0.11803456, 0.37568667 }, 
	{ 0.47357251, 0.11501071, 0.38338399 }, 
	{ 0.47924506, 0.11186553, 0.39130864 }, 
	{ 0.48482535, 0.10861132, 0.39945883 }, 
	{ 0.49030394, 0.10526086, 0.40783299 }, 
	{ 0.49566975, 0.10183585, 0.4164232  }, 
	{ 0.50091274, 0.09835974, 0.4252212  }, 
	{ 0.50602315, 0.09486118, 0.43421594 }, 
	{ 0.51099169, 0.0913744 , 0.4433937  }, 
	{ 0.51580979, 0.0879393 , 0.45273836 }, 
	{ 0.52046983, 0.08460151, 0.46223172 }, 
	{ 0.52496528, 0.08141217, 0.47185394 }, 
	{ 0.52929077, 0.07842746, 0.48158397 }, 
	{ 0.53344214, 0.07570762, 0.49140009 }, 
	{ 0.53741646, 0.07331545, 0.50128038 }, 
	{ 0.54121191, 0.07131449, 0.51120295 }, 
	{ 0.54482771, 0.0697663 , 0.52114633 }, 
	{ 0.54826405, 0.068727  , 0.53109007 }, 
	{ 0.5515219 , 0.06824411, 0.54101477 }, 
	{ 0.55460291, 0.06835358, 0.55090216 }, 
	{ 0.55750924, 0.06907723, 0.56073534 }, 
	{ 0.56024349, 0.07042163, 0.57049884 }, 
	{ 0.56280851, 0.07237845, 0.58017846 }, 
	{ 0.5652074 , 0.074926  , 0.58976132 }, 
	{ 0.56744336, 0.07803182, 0.59923572 }, 
	{ 0.56951969, 0.08165576, 0.60859094 }, 
	{ 0.57143963, 0.08575302, 0.61781771 }, 
	{ 0.57320642, 0.09027707, 0.62690754 }, 
	{ 0.57482327, 0.09518184, 0.63585281 }, 
	{ 0.57629328, 0.1004233 , 0.64464669 }, 
	{ 0.57761947, 0.1059605 , 0.65328307 }, 
	{ 0.57880479, 0.11175607, 0.66175646 }, 
	{ 0.57985206, 0.11777648, 0.67006198 }, 
	{ 0.58076404, 0.1239919 , 0.67819522 }, 
	{ 0.58154339, 0.13037605, 0.68615227 }, 
	{ 0.58219268, 0.13690587, 0.6939296  }, 
	{ 0.58271442, 0.1435612 , 0.70152406 }, 
	{ 0.58311105, 0.15032442, 0.70893284 }, 
	{ 0.58338497, 0.15718015, 0.71615342 }, 
	{ 0.58353852, 0.16411493, 0.72318358 }, 
	{ 0.58357402, 0.17111698, 0.73002132 }, 
	{ 0.58349377, 0.17817593, 0.7366649  }, 
	{ 0.58330008, 0.18528263, 0.74311279 }, 
	{ 0.58299524, 0.19242899, 0.74936366 }, 
	{ 0.58258158, 0.19960777, 0.75541636 }, 
	{ 0.58206142, 0.20681253, 0.76126998 }, 
	{ 0.58143712, 0.21403747, 0.76692374 }, 
	{ 0.58071112, 0.22127732, 0.77237705 }, 
	{ 0.57988591, 0.22852728, 0.77762948 }, 
	{ 0.57896402, 0.23578296, 0.78268078 }, 
	{ 0.57794809, 0.2430403 , 0.78753087 }, 
	{ 0.57684084, 0.25029554, 0.79217983 }, 
	{ 0.57564508, 0.25754517, 0.7966279  }, 
	{ 0.57436373, 0.26478589, 0.80087551 }, 
	{ 0.57299981, 0.27201459, 0.80492325 }, 
	{ 0.57155649, 0.27922832, 0.80877191 }, 
	{ 0.57003705, 0.28642427, 0.81242242 }, 
	{ 0.56844491, 0.29359975, 0.81587594 }, 
	{ 0.56678362, 0.30075218, 0.81913378 }, 
	{ 0.56505689, 0.30787908, 0.82219746 }, 
	{ 0.56326858, 0.31497804, 0.8250687  }, 
	{ 0.56142269, 0.32204677, 0.8277494  }, 
	{ 0.55952336, 0.32908303, 0.8302417  }, 
	{ 0.55757491, 0.33608464, 0.83254789 }, 
	{ 0.55558182, 0.34304953, 0.83467053 }, 
	{ 0.55354871, 0.34997564, 0.83661233 }, 
	{ 0.55148034, 0.35686103, 0.83837625 }, 
	{ 0.54938162, 0.36370381, 0.83996545 }, 
	{ 0.54725759, 0.37050216, 0.84138328 }, 
	{ 0.54511344, 0.37725433, 0.84263332 }, 
	{ 0.54295446, 0.38395866, 0.84371932 }, 
	{ 0.54078604, 0.39061356, 0.84464526 }, 
	{ 0.53861369, 0.39721751, 0.84541529 }, 
	{ 0.53644299, 0.4037691 , 0.84603374 }, 
	{ 0.53427957, 0.41026698, 0.84650513 }, 
	{ 0.53212915, 0.41670992, 0.84683411 }, 
	{ 0.52999745, 0.42309676, 0.84702553 }, 
	{ 0.52789021, 0.42942645, 0.84708434 }, 
	{ 0.52581315, 0.43569805, 0.84701564 }, 
	{ 0.523772  , 0.4419107 , 0.84682463 }, 
	{ 0.5217724 , 0.44806368, 0.84651661 }, 
	{ 0.51981993, 0.45415634, 0.84609697 }, 
	{ 0.51792009, 0.46018816, 0.84557116 }, 
	{ 0.51607826, 0.46615872, 0.84494467 }, 
	{ 0.51429966, 0.47206772, 0.84422304 }, 
	{ 0.5125894 , 0.47791495, 0.84341181 }, 
	{ 0.51095236, 0.48370033, 0.84251655 }, 
	{ 0.50939329, 0.48942384, 0.84154281 }, 
	{ 0.50791666, 0.49508561, 0.84049609 }, 
	{ 0.50652672, 0.50068587, 0.83938183 }, 
	{ 0.50522749, 0.50622492, 0.83820544 }, 
	{ 0.50402271, 0.51170318, 0.83697223 }, 
	{ 0.50291585, 0.51712114, 0.83568745 }, 
	{ 0.5019101 , 0.52247939, 0.83435622 }, 
	{ 0.50100834, 0.52777859, 0.83298359 }, 
	{ 0.50021319, 0.53301947, 0.83157449 }, 
	{ 0.49952693, 0.53820284, 0.83013373 }, 
	{ 0.4989515 , 0.5433296 , 0.8286659  }, 
	{ 0.49848855, 0.54840071, 0.82717552 }, 
	{ 0.49813941, 0.55341717, 0.82566693 }, 
	{ 0.49790512, 0.55838004, 0.82414434 }, 
	{ 0.49778644, 0.5632904 , 0.82261188 }, 
	{ 0.49778378, 0.5681494 , 0.82107343 }, 
	{ 0.49789723, 0.57295826, 0.81953267 }, 
	{ 0.49812663, 0.57771821, 0.81799316 }, 
	{ 0.49847153, 0.58243051, 0.81645831 }, 
	{ 0.49893128, 0.58709639, 0.8149315  }, 
	{ 0.49950487, 0.59171721, 0.81341573 }, 
	{ 0.50019106, 0.59629429, 0.81191388 }, 
	{ 0.50098841, 0.600829  , 0.81042874 }, 
	{ 0.50189529, 0.60532265, 0.80896303 }, 
	{ 0.50290981, 0.60977665, 0.80751914 }, 
	{ 0.50402985, 0.61419238, 0.80609933 }, 
	{ 0.50525318, 0.61857124, 0.80470576 }, 
	{ 0.50657743, 0.62291459, 0.80334059 }, 
	{ 0.50799996, 0.62722387, 0.80200552 }, 
	{ 0.50951805, 0.63150048, 0.80070228 }, 
	{ 0.51112889, 0.63574581, 0.79943254 }, 
	{ 0.51282948, 0.63996127, 0.79819768 }, 
	{ 0.51461675, 0.64414827, 0.79699896 }, 
	{ 0.51648755, 0.64830822, 0.79583758 }, 
	{ 0.51843862, 0.6524425 , 0.79471458 }, 
	{ 0.52046664, 0.65655254, 0.79363078 }, 
	{ 0.52256822, 0.66063972, 0.79258695 }, 
	{ 0.52473993, 0.66470545, 0.79158372 }, 
	{ 0.52697828, 0.66875111, 0.79062152 }, 
	{ 0.52927978, 0.67277809, 0.78970067 }, 
	{ 0.53164088, 0.67678779, 0.78882137 }, 
	{ 0.53405803, 0.68078158, 0.78798362 }, 
	{ 0.53652767, 0.68476085, 0.78718733 }, 
	{ 0.53904627, 0.68872697, 0.7864322  }, 
	{ 0.54161027, 0.69268131, 0.78571781 }, 
	{ 0.54421616, 0.69662523, 0.78504357 }, 
	{ 0.54686046, 0.7005601 , 0.78440868 }, 
	{ 0.54953973, 0.70448725, 0.78381224 }, 
	{ 0.5522506 , 0.70840802, 0.78325313 }, 
	{ 0.55498973, 0.71232375, 0.78272997 }, 
	{ 0.55775391, 0.71623573, 0.78224135 }, 
	{ 0.56054   , 0.72014527, 0.78178559 }, 
	{ 0.56334495, 0.72405364, 0.78136067 }, 
	{ 0.56616586, 0.72796208, 0.7809646  }, 
	{ 0.56899999, 0.73187181, 0.78059509 }, 
	{ 0.57184471, 0.73578401, 0.78024954 }, 
	{ 0.5746976 , 0.73969985, 0.77992519 }, 
	{ 0.57755644, 0.7436204 , 0.77961919 }, 
	{ 0.58041922, 0.74754671, 0.77932838 }, 
	{ 0.58328416, 0.75147981, 0.77904926 }, 
	{ 0.58614975, 0.75542061, 0.7787783  }, 
	{ 0.58901475, 0.75936996, 0.7785118  }, 
	{ 0.59187822, 0.76332865, 0.7782458  }, 
	{ 0.59473954, 0.76729738, 0.77797618 }, 
	{ 0.59759841, 0.77127679, 0.7776985  }, 
	{ 0.60045488, 0.77526737, 0.77740843 }, 
	{ 0.60330938, 0.77926953, 0.77710145 }, 
	{ 0.60616272, 0.78328359, 0.77677292 }, 
	{ 0.60901606, 0.78730972, 0.77641811 }, 
	{ 0.61187097, 0.79134802, 0.77603224 }, 
	{ 0.61472943, 0.79539844, 0.77561052 }, 
	{ 0.61759377, 0.79946082, 0.77514814 }, 
	{ 0.62046677, 0.80353487, 0.77464034 }, 
	{ 0.62335154, 0.80762017, 0.77408242 }, 
	{ 0.62625161, 0.8117162 , 0.77346977 }, 
	{ 0.62917087, 0.81582228, 0.7727979  }, 
	{ 0.63211355, 0.81993764, 0.77206245 }, 
	{ 0.63508424, 0.82406137, 0.77125927 }, 
	{ 0.63808784, 0.82819246, 0.77038427 }, 
	{ 0.64112954, 0.83232981, 0.76943331 }, 
	{ 0.64421489, 0.83647212, 0.76840325 }, 
	{ 0.64734968, 0.84061803, 0.76729084 }, 
	{ 0.65053992, 0.84476616, 0.76609234 }, 
	{ 0.65379192, 0.84891487, 0.76480547 }, 
	{ 0.6571122 , 0.85306252, 0.76342757 }, 
	{ 0.66050751, 0.85720734, 0.76195622 }, 
	{ 0.6639848 , 0.86134747, 0.7603896  }, 
	{ 0.66755124, 0.86548088, 0.75872646 }, 
	{ 0.67121422, 0.86960558, 0.75696447 }, 
	{ 0.67498129, 0.87371931, 0.75510363 }, 
	{ 0.67886025, 0.87781976, 0.75314325 }, 
	{ 0.68285911, 0.8819045 , 0.75108309 }, 
	{ 0.6869861 , 0.88597096, 0.7489234  }, 
	{ 0.69124973, 0.89001641, 0.74666498 }, 
	{ 0.69565872, 0.89403798, 0.74430926 }, 
	{ 0.70022207, 0.89803257, 0.74185843 }, 
	{ 0.70494935, 0.90199697, 0.73931381 }, 
	{ 0.70984994, 0.90592765, 0.73668045 }, 
	{ 0.71493422, 0.90982091, 0.73396099 }, 
	{ 0.72021236, 0.91367274, 0.73116189 }, 
	{ 0.72569528, 0.91747884, 0.72828927 }, 
	{ 0.7313942 , 0.9212346 , 0.72535097 }, 
	{ 0.73732063, 0.92493506, 0.72235684 }, 
	{ 0.74348624, 0.92857492, 0.71931896 }, 
	{ 0.74990347, 0.93214846, 0.71625026 }, 
	{ 0.7565837 , 0.93564971, 0.71316899 }, 
	{ 0.76353852, 0.93907236, 0.71009516 }, 
	{ 0.77077848, 0.94240996, 0.70705291 }, 
	{ 0.77831226, 0.94565608, 0.70407107 }, 
	{ 0.78614703, 0.9488044 , 0.70118136 }, 
	{ 0.79428492, 0.9518495 , 0.69842228 }, 
	{ 0.80272416, 0.95478691, 0.69583468 }, 
	{ 0.81145663, 0.9576139 , 0.69346237 }, 
	{ 0.82046728, 0.96032991, 0.69134936 }, 
	{ 0.82973371, 0.96293708, 0.68953698 }, 
	{ 0.83922653, 0.9654405 , 0.68806037 }, 
	{ 0.8489114 , 0.96784792, 0.68694452 }, 
	{ 0.85875138, 0.97016928, 0.68620163 }, 
	{ 0.86870981, 0.97241593, 0.68582981 }, 
	{ 0.87875336, 0.97459951, 0.68581323 }
};

// contributed by @cmj
float inferno[256][3] = {
	{ 0.001462, 0.000466, 0.013866 },
	{ 0.002267, 0.001270, 0.018570 },
	{ 0.003299, 0.002249, 0.024239 },
	{ 0.004547, 0.003392, 0.030909 },
	{ 0.006006, 0.004692, 0.038558 },
	{ 0.007676, 0.006136, 0.046836 },
	{ 0.009561, 0.007713, 0.055143 },
	{ 0.011663, 0.009417, 0.063460 },
	{ 0.013995, 0.011225, 0.071862 },
	{ 0.016561, 0.013136, 0.080282 },
	{ 0.019373, 0.015133, 0.088767 },
	{ 0.022447, 0.017199, 0.097327 },
	{ 0.025793, 0.019331, 0.105930 },
	{ 0.029432, 0.021503, 0.114621 },
	{ 0.033385, 0.023702, 0.123397 },
	{ 0.037668, 0.025921, 0.132232 },
	{ 0.042253, 0.028139, 0.141141 },
	{ 0.046915, 0.030324, 0.150164 },
	{ 0.051644, 0.032474, 0.159254 },
	{ 0.056449, 0.034569, 0.168414 },
	{ 0.061340, 0.036590, 0.177642 },
	{ 0.066331, 0.038504, 0.186962 },
	{ 0.071429, 0.040294, 0.196354 },
	{ 0.076637, 0.041905, 0.205799 },
	{ 0.081962, 0.043328, 0.215289 },
	{ 0.087411, 0.044556, 0.224813 },
	{ 0.092990, 0.045583, 0.234358 },
	{ 0.098702, 0.046402, 0.243904 },
	{ 0.104551, 0.047008, 0.253430 },
	{ 0.110536, 0.047399, 0.262912 },
	{ 0.116656, 0.047574, 0.272321 },
	{ 0.122908, 0.047536, 0.281624 },
	{ 0.129285, 0.047293, 0.290788 },
	{ 0.135778, 0.046856, 0.299776 },
	{ 0.142378, 0.046242, 0.308553 },
	{ 0.149073, 0.045468, 0.317085 },
	{ 0.155850, 0.044559, 0.325338 },
	{ 0.162689, 0.043554, 0.333277 },
	{ 0.169575, 0.042489, 0.340874 },
	{ 0.176493, 0.041402, 0.348111 },
	{ 0.183429, 0.040329, 0.354971 },
	{ 0.190367, 0.039309, 0.361447 },
	{ 0.197297, 0.038400, 0.367535 },
	{ 0.204209, 0.037632, 0.373238 },
	{ 0.211095, 0.037030, 0.378563 },
	{ 0.217949, 0.036615, 0.383522 },
	{ 0.224763, 0.036405, 0.388129 },
	{ 0.231538, 0.036405, 0.392400 },
	{ 0.238273, 0.036621, 0.396353 },
	{ 0.244967, 0.037055, 0.400007 },
	{ 0.251620, 0.037705, 0.403378 },
	{ 0.258234, 0.038571, 0.406485 },
	{ 0.264810, 0.039647, 0.409345 },
	{ 0.271347, 0.040922, 0.411976 },
	{ 0.277850, 0.042353, 0.414392 },
	{ 0.284321, 0.043933, 0.416608 },
	{ 0.290763, 0.045644, 0.418637 },
	{ 0.297178, 0.047470, 0.420491 },
	{ 0.303568, 0.049396, 0.422182 },
	{ 0.309935, 0.051407, 0.423721 },
	{ 0.316282, 0.053490, 0.425116 },
	{ 0.322610, 0.055634, 0.426377 },
	{ 0.328921, 0.057827, 0.427511 },
	{ 0.335217, 0.060060, 0.428524 },
	{ 0.341500, 0.062325, 0.429425 },
	{ 0.347771, 0.064616, 0.430217 },
	{ 0.354032, 0.066925, 0.430906 },
	{ 0.360284, 0.069247, 0.431497 },
	{ 0.366529, 0.071579, 0.431994 },
	{ 0.372768, 0.073915, 0.432400 },
	{ 0.379001, 0.076253, 0.432719 },
	{ 0.385228, 0.078591, 0.432955 },
	{ 0.391453, 0.080927, 0.433109 },
	{ 0.397674, 0.083257, 0.433183 },
	{ 0.403894, 0.085580, 0.433179 },
	{ 0.410113, 0.087896, 0.433098 },
	{ 0.416331, 0.090203, 0.432943 },
	{ 0.422549, 0.092501, 0.432714 },
	{ 0.428768, 0.094790, 0.432412 },
	{ 0.434987, 0.097069, 0.432039 },
	{ 0.441207, 0.099338, 0.431594 },
	{ 0.447428, 0.101597, 0.431080 },
	{ 0.453651, 0.103848, 0.430498 },
	{ 0.459875, 0.106089, 0.429846 },
	{ 0.466100, 0.108322, 0.429125 },
	{ 0.472328, 0.110547, 0.428334 },
	{ 0.478558, 0.112764, 0.427475 },
	{ 0.484789, 0.114974, 0.426548 },
	{ 0.491022, 0.117179, 0.425552 },
	{ 0.497257, 0.119379, 0.424488 },
	{ 0.503493, 0.121575, 0.423356 },
	{ 0.509730, 0.123769, 0.422156 },
	{ 0.515967, 0.125960, 0.420887 },
	{ 0.522206, 0.128150, 0.419549 },
	{ 0.528444, 0.130341, 0.418142 },
	{ 0.534683, 0.132534, 0.416667 },
	{ 0.540920, 0.134729, 0.415123 },
	{ 0.547157, 0.136929, 0.413511 },
	{ 0.553392, 0.139134, 0.411829 },
	{ 0.559624, 0.141346, 0.410078 },
	{ 0.565854, 0.143567, 0.408258 },
	{ 0.572081, 0.145797, 0.406369 },
	{ 0.578304, 0.148039, 0.404411 },
	{ 0.584521, 0.150294, 0.402385 },
	{ 0.590734, 0.152563, 0.400290 },
	{ 0.596940, 0.154848, 0.398125 },
	{ 0.603139, 0.157151, 0.395891 },
	{ 0.609330, 0.159474, 0.393589 },
	{ 0.615513, 0.161817, 0.391219 },
	{ 0.621685, 0.164184, 0.388781 },
	{ 0.627847, 0.166575, 0.386276 },
	{ 0.633998, 0.168992, 0.383704 },
	{ 0.640135, 0.171438, 0.381065 },
	{ 0.646260, 0.173914, 0.378359 },
	{ 0.652369, 0.176421, 0.375586 },
	{ 0.658463, 0.178962, 0.372748 },
	{ 0.664540, 0.181539, 0.369846 },
	{ 0.670599, 0.184153, 0.366879 },
	{ 0.676638, 0.186807, 0.363849 },
	{ 0.682656, 0.189501, 0.360757 },
	{ 0.688653, 0.192239, 0.357603 },
	{ 0.694627, 0.195021, 0.354388 },
	{ 0.700576, 0.197851, 0.351113 },
	{ 0.706500, 0.200728, 0.347777 },
	{ 0.712396, 0.203656, 0.344383 },
	{ 0.718264, 0.206636, 0.340931 },
	{ 0.724103, 0.209670, 0.337424 },
	{ 0.729909, 0.212759, 0.333861 },
	{ 0.735683, 0.215906, 0.330245 },
	{ 0.741423, 0.219112, 0.326576 },
	{ 0.747127, 0.222378, 0.322856 },
	{ 0.752794, 0.225706, 0.319085 },
	{ 0.758422, 0.229097, 0.315266 },
	{ 0.764010, 0.232554, 0.311399 },
	{ 0.769556, 0.236077, 0.307485 },
	{ 0.775059, 0.239667, 0.303526 },
	{ 0.780517, 0.243327, 0.299523 },
	{ 0.785929, 0.247056, 0.295477 },
	{ 0.791293, 0.250856, 0.291390 },
	{ 0.796607, 0.254728, 0.287264 },
	{ 0.801871, 0.258674, 0.283099 },
	{ 0.807082, 0.262692, 0.278898 },
	{ 0.812239, 0.266786, 0.274661 },
	{ 0.817341, 0.270954, 0.270390 },
	{ 0.822386, 0.275197, 0.266085 },
	{ 0.827372, 0.279517, 0.261750 },
	{ 0.832299, 0.283913, 0.257383 },
	{ 0.837165, 0.288385, 0.252988 },
	{ 0.841969, 0.292933, 0.248564 },
	{ 0.846709, 0.297559, 0.244113 },
	{ 0.851384, 0.302260, 0.239636 },
	{ 0.855992, 0.307038, 0.235133 },
	{ 0.860533, 0.311892, 0.230606 },
	{ 0.865006, 0.316822, 0.226055 },
	{ 0.869409, 0.321827, 0.221482 },
	{ 0.873741, 0.326906, 0.216886 },
	{ 0.878001, 0.332060, 0.212268 },
	{ 0.882188, 0.337287, 0.207628 },
	{ 0.886302, 0.342586, 0.202968 },
	{ 0.890341, 0.347957, 0.198286 },
	{ 0.894305, 0.353399, 0.193584 },
	{ 0.898192, 0.358911, 0.188860 },
	{ 0.902003, 0.364492, 0.184116 },
	{ 0.905735, 0.370140, 0.179350 },
	{ 0.909390, 0.375856, 0.174563 },
	{ 0.912966, 0.381636, 0.169755 },
	{ 0.916462, 0.387481, 0.164924 },
	{ 0.919879, 0.393389, 0.160070 },
	{ 0.923215, 0.399359, 0.155193 },
	{ 0.926470, 0.405389, 0.150292 },
	{ 0.929644, 0.411479, 0.145367 },
	{ 0.932737, 0.417627, 0.140417 },
	{ 0.935747, 0.423831, 0.135440 },
	{ 0.938675, 0.430091, 0.130438 },
	{ 0.941521, 0.436405, 0.125409 },
	{ 0.944285, 0.442772, 0.120354 },
	{ 0.946965, 0.449191, 0.115272 },
	{ 0.949562, 0.455660, 0.110164 },
	{ 0.952075, 0.462178, 0.105031 },
	{ 0.954506, 0.468744, 0.099874 },
	{ 0.956852, 0.475356, 0.094695 },
	{ 0.959114, 0.482014, 0.089499 },
	{ 0.961293, 0.488716, 0.084289 },
	{ 0.963387, 0.495462, 0.079073 },
	{ 0.965397, 0.502249, 0.073859 },
	{ 0.967322, 0.509078, 0.068659 },
	{ 0.969163, 0.515946, 0.063488 },
	{ 0.970919, 0.522853, 0.058367 },
	{ 0.972590, 0.529798, 0.053324 },
	{ 0.974176, 0.536780, 0.048392 },
	{ 0.975677, 0.543798, 0.043618 },
	{ 0.977092, 0.550850, 0.039050 },
	{ 0.978422, 0.557937, 0.034931 },
	{ 0.979666, 0.565057, 0.031409 },
	{ 0.980824, 0.572209, 0.028508 },
	{ 0.981895, 0.579392, 0.026250 },
	{ 0.982881, 0.586606, 0.024661 },
	{ 0.983779, 0.593849, 0.023770 },
	{ 0.984591, 0.601122, 0.023606 },
	{ 0.985315, 0.608422, 0.024202 },
	{ 0.985952, 0.615750, 0.025592 },
	{ 0.986502, 0.623105, 0.027814 },
	{ 0.986964, 0.630485, 0.030908 },
	{ 0.987337, 0.637890, 0.034916 },
	{ 0.987622, 0.645320, 0.039886 },
	{ 0.987819, 0.652773, 0.045581 },
	{ 0.987926, 0.660250, 0.051750 },
	{ 0.987945, 0.667748, 0.058329 },
	{ 0.987874, 0.675267, 0.065257 },
	{ 0.987714, 0.682807, 0.072489 },
	{ 0.987464, 0.690366, 0.079990 },
	{ 0.987124, 0.697944, 0.087731 },
	{ 0.986694, 0.705540, 0.095694 },
	{ 0.986175, 0.713153, 0.103863 },
	{ 0.985566, 0.720782, 0.112229 },
	{ 0.984865, 0.728427, 0.120785 },
	{ 0.984075, 0.736087, 0.129527 },
	{ 0.983196, 0.743758, 0.138453 },
	{ 0.982228, 0.751442, 0.147565 },
	{ 0.981173, 0.759135, 0.156863 },
	{ 0.980032, 0.766837, 0.166353 },
	{ 0.978806, 0.774545, 0.176037 },
	{ 0.977497, 0.782258, 0.185923 },
	{ 0.976108, 0.789974, 0.196018 },
	{ 0.974638, 0.797692, 0.206332 },
	{ 0.973088, 0.805409, 0.216877 },
	{ 0.971468, 0.813122, 0.227658 },
	{ 0.969783, 0.820825, 0.238686 },
	{ 0.968041, 0.828515, 0.249972 },
	{ 0.966243, 0.836191, 0.261534 },
	{ 0.964394, 0.843848, 0.273391 },
	{ 0.962517, 0.851476, 0.285546 },
	{ 0.960626, 0.859069, 0.298010 },
	{ 0.958720, 0.866624, 0.310820 },
	{ 0.956834, 0.874129, 0.323974 },
	{ 0.954997, 0.881569, 0.337475 },
	{ 0.953215, 0.888942, 0.351369 },
	{ 0.951546, 0.896226, 0.365627 },
	{ 0.950018, 0.903409, 0.380271 },
	{ 0.948683, 0.910473, 0.395289 },
	{ 0.947594, 0.917399, 0.410665 },
	{ 0.946809, 0.924168, 0.426373 },
	{ 0.946392, 0.930761, 0.442367 },
	{ 0.946403, 0.937159, 0.458592 },
	{ 0.946903, 0.943348, 0.474970 },
	{ 0.947937, 0.949318, 0.491426 },
	{ 0.949545, 0.955063, 0.507860 },
	{ 0.951740, 0.960587, 0.524203 },
	{ 0.954529, 0.965896, 0.540361 },
	{ 0.957896, 0.971003, 0.556275 },
	{ 0.961812, 0.975924, 0.571925 },
	{ 0.966249, 0.980678, 0.587206 },
	{ 0.971162, 0.985282, 0.602154 },
	{ 0.976511, 0.989753, 0.616760 },
	{ 0.982257, 0.994109, 0.631017 },
	{ 0.988362, 0.998364, 0.644924 }
};

void CPlotter::setWfColormap(const QString &cmap)
{
    int i;

    if (cmap.compare("gqrx", Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
        {
            // level 0: black background
            if (i < 20)
                m_ColorTbl[i].setRgb(0, 0, 0);
                // level 1: black -> blue
            else if ((i >= 20) && (i < 70))
                m_ColorTbl[i].setRgb(0, 0, 140*(i-20)/50);
                // level 2: blue -> light-blue / greenish
            else if ((i >= 70) && (i < 100))
                m_ColorTbl[i].setRgb(60*(i-70)/30, 125*(i-70)/30, 115*(i-70)/30 + 140);
                // level 3: light blue -> yellow
            else if ((i >= 100) && (i < 150))
                m_ColorTbl[i].setRgb(195*(i-100)/50 + 60, 130*(i-100)/50 + 125, 255-(255*(i-100)/50));
                // level 4: yellow -> red
            else if ((i >= 150) && (i < 250))
                m_ColorTbl[i].setRgb(255, 255-255*(i-150)/100, 0);
                // level 5: red -> white
            else if (i >= 250)
                m_ColorTbl[i].setRgb(255, 255*(i-250)/5, 255*(i-250)/5);
        }
    }
    else if (cmap.compare("turbo", Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(turbo[i][0], turbo[i][1], turbo[i][2]);
    }
    else if (cmap.compare("plasma",Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(plasma[i][0], plasma[i][1], plasma[i][2]);
    }
    else if (cmap.compare("whitehotcompressed",Qt::CaseInsensitive) == 0)
    {
        // contributed by @drmpeg @devnulling
        // for use with high quality spectrum paining
        // see https://gist.github.com/drmpeg/31a9a7dd6918856aeb60
        for (i = 0; i < 256; i++)
        {
            if (i < 64)
            {
                m_ColorTbl[i].setRgb(i*4, i*4, i*4);
            }
            else
            {
                m_ColorTbl[i].setRgb(255, 255, 255);
            }
        }
    }
    else if (cmap.compare("whitehot",Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(i, i, i);
    }
    else if (cmap.compare("blackhot",Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(255-i, 255-i, 255-i);
    }
    else if (cmap.compare("viridis",Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(viridis[i][0] * 256, viridis[i][1] * 256, viridis[i][2] * 256);
    }
    else if (cmap.compare("magma",Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(magma[i][0] * 256, magma[i][1] * 256, magma[i][2] * 256);
    }
    else if (cmap.compare("inferno",Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(inferno[i][0] * 256, inferno[i][1] * 256, inferno[i][2] * 256);
    }
    else if (cmap.compare("grape",Qt::CaseInsensitive) == 0)
    {
        for (i = 0; i < 256; i++)
            m_ColorTbl[i].setRgb(grape[i][0] * 256, grape[i][1] * 256, grape[i][2] * 256);
    }

}
