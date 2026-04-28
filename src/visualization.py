import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


def animate_simulation(
    history,
    interval_ms=30,
    cmap="RdBu_r",
    lens_mask=None,
    save_path=None,
    show=True,
    fps=None,
    dpi=120,
):
    if len(history) == 0:
        raise ValueError("history is empty; run the simulation with steps > 0")

    fig, ax = plt.subplots()
    use_blit = True

    max_abs = max(float(np.max(np.abs(frame))) for frame in history)
    if max_abs == 0.0:
        max_abs = 1.0

    im = ax.imshow(
        history[0].T,
        cmap=cmap,
        origin="lower",
        vmin=-max_abs,
        vmax=max_abs,
        animated=True,
    )
    title = ax.set_title("2D Wave Equation (Finite Differences)")

    lens_artists = []
    if lens_mask is not None:
        lens_mask = np.asarray(lens_mask)
        if lens_mask.shape != history[0].shape:
            raise ValueError(
                f"lens_mask shape {lens_mask.shape} must match frame shape {history[0].shape}"
            )

        # Frames are plotted as history[...].T, so transpose the mask as well.
        cs = ax.contour(
            lens_mask.T.astype(float),
            levels=[0.5],
            colors="k",
            linewidths=1.5,
            alpha=0.9,
        )
        # Matplotlib API difference across versions:
        # some expose the drawn artists as `collections`, others don't.
        if hasattr(cs, "collections"):
            lens_artists = list(cs.collections)
            for coll in lens_artists:
                coll.set_animated(True)
        else:
            # Fallback: keep the contour static and disable blitting so it stays visible.
            use_blit = False

    def update(frame):
        im.set_array(history[frame].T)
        title.set_text(f"2D Wave Equation (Frame {frame + 1}/{len(history)})")
        if use_blit:
            return [im, title, *lens_artists]
        return [im, title]

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(history),
        interval=interval_ms,
        blit=use_blit,
    )

    plt.colorbar(im, ax=ax)

    if fps is None:
        fps = 1000.0 / float(interval_ms)

    if save_path is not None:
        save_path = str(save_path)
        lower = save_path.lower()
        if lower.endswith(".gif"):
            writer = animation.PillowWriter(fps=fps)
        elif lower.endswith(".mp4"):
            writer = animation.FFMpegWriter(fps=fps)
        else:
            # Let Matplotlib pick the best available writer for unknown extensions.
            writer = None

        if writer is None:
            ani.save(save_path, dpi=dpi)
        else:
            ani.save(save_path, writer=writer, dpi=dpi)

    if show:
        plt.show()
    else:
        plt.close(fig)
