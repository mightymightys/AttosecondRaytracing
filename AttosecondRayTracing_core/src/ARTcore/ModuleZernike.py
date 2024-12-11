from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


def zernike_gradient(x, y, max_order):
    """Recursively generate Zernike and Zernike Gradient polynomials and evaluate at input x and y coordinates
        Algorithm and background from this paper:
        Torben B. Andersen, "Efficient and robust recurrence relations for the Zernike circle polynomials and their derivatives in Cartesian coordinates,"
        Opt. Express 26, 18878-18896 (2018)

    // Author : Logan Rodriguez Graves. Date: 1/27/022


        Parameters
        ----------
        x : numpy.ndarray
            x coordinates to evaluate at
        y : numpy.ndarray
            y coordinates to evaluate at
        max_order : int
            polynomial order

        Returns
        -------
        Note: n, m order is the radial and azimuthal order respectively, with m = 0:n
        dict
            zernike polynomial evaluated at the given points, with key values of the n,m order

        dict
            zernike x gradient  polynomial evaluated at the given points, with key values of the n, m order

        dict
            zernike y gradient polynomial evaluated at the given points, with key values of the n, m order
    """

    # Determine the max radial order, if it is less than 2 then set it to at least 2.
    if max_order < 2:
        max_order = 2

    num_coords = len(x)

    # Dict is your friend here! They will allow us to directly reference into
    # the sets we need, as opposed to an unnecessarily extra complex step of translating nm
    # indices to a single order.

    # declare a default dict
    zernike_map = defaultdict(list)
    gradient_x_map = defaultdict(list)
    gradient_y_map = defaultdict(list)

    # Assign seed values
    zernike_map[(0, 0)].append((0, [1 for i in range(len(x))]))
    zernike_map[(1, 0)].append((1, y))
    zernike_map[(1, 1)].append((2, x))

    # Assign seed gradient values
    gradient_x_map[(0, 0)].append((0, [0 for i in range(len(x))]))
    gradient_x_map[(1, 0)].append((1, [0 for i in range(len(x))]))
    gradient_x_map[(1, 1)].append((2, [1 for i in range(len(x))]))

    gradient_y_map[(0, 0)].append((0, [0 for i in range(len(x))]))
    gradient_y_map[(1, 0)].append((1, [1 for i in range(len(x))]))
    gradient_y_map[(1, 1)].append((2, [0 for i in range(len(x))]))

    # Define a variable for the ansi order going forward
    ansi_order = 2

    # Outer loop for radial order index
    for n in range(2, max_order + 1):
        # Inner loop for azimuthal index
        for m in range(0, n + 1):
            # Step the ansi order for our next iterations
            ansi_order += 1

            new_zern_values = []
            new_gradient_x = []
            new_gradient_y = []

            # Handle the various cases of exception
            if m == 0:
                # Collect Needed Prior Zernike and Gradient Polynomials
                # Index(n - 1, 0)
                zern_n1_0 = zernike_map.get((n - 1, 0))[0][1]
                # Index(n - 1, n - 1)
                zern_n1_n1 = zernike_map.get((n - 1, n - 1))[0][1]

                # Iterating through every input coordinate point, define the new
                # Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in range(0, num_coords):
                    new_zern_values.append(x[pos] * zern_n1_0[pos] + y[pos] * zern_n1_n1[pos])

                    new_gradient_x.append(n * zern_n1_0[pos])

                    new_gradient_y.append(n * zern_n1_n1[pos])

            elif m == n:
                # Collect Needed Prior Zernike and Gradient Polynomials
                # Index(n - 1, 0)
                zern_n1_0 = zernike_map.get((n - 1, 0))[0][1]
                # Index(n - 1, n - 1)
                zern_n1_n1 = zernike_map.get((n - 1, n - 1))[0][1]

                # Iterating through every input coordinate point, define the new
                # Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in range(0, num_coords):
                    new_zern_values.append(x[pos] * zern_n1_n1[pos] - y[pos] * zern_n1_0[pos])
                    new_gradient_x.append(n * zern_n1_n1[pos])
                    new_gradient_y.append(-1.0 * n * zern_n1_0[pos])

            elif n % 2 != 0 and m == (n - 1) / 2:
                # Collect Needed Prior Zernike and Gradient Polynomials
                # Index(n - 1, n - 1 - m)
                zern_n1_n1m = zernike_map.get((n - 1, n - 1 - m))[0][1]

                # Index(n - 1, m - 1)
                zern_n1_m1 = zernike_map.get((n - 1, m - 1))[0][1]

                # Index(n - 1, n - m)
                zern_n1_nm = zernike_map.get((n - 1, n - m))[0][1]

                # Index(n - 2, m - 1)
                zern_n2_m1 = zernike_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_x_n2_m1 = gradient_x_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_y_n2_m1 = gradient_y_map.get((n - 2, m - 1))[0][1]

                # Iterating through every input coordinate point, define the new
                # Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in range(0, num_coords):
                    # Calculate new values
                    new_zern_values.append(
                        y[pos] * zern_n1_n1m[pos]
                        + x[pos] * zern_n1_m1[pos]
                        - y[pos] * zern_n1_nm[pos]
                        - zern_n2_m1[pos]
                    )

                    new_gradient_x.append(n * zern_n1_m1[pos] + gradient_x_n2_m1[pos])
                    new_gradient_y.append(n * zern_n1_n1m[pos] - n * zern_n1_nm[pos] + gradient_y_n2_m1[pos])

            elif n % 2 != 0 and m == (n - 1) / 2 + 1:
                # Collect Needed Prior Zernike and Gradient Polynomials

                # Index(n - 1, n - 1 - m)
                zern_n1_n1m = zernike_map.get((n - 1, n - 1 - m))[0][1]

                # Index(n - 1, m)
                zern_n1_m = zernike_map.get((n - 1, m))[0][1]

                # Index(n - 1, m - 1)
                zern_n1_m1 = zernike_map.get((n - 1, m - 1))[0][1]

                # Index(n - 2, m - 1)
                zern_n2_m1 = zernike_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_x_n2_m1 = gradient_x_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_y_n2_m1 = gradient_y_map.get((n - 2, m - 1))[0][1]

                # Iterating through every input coordinate point, define the new
                # Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in range(0, num_coords):
                    # Calculate new values
                    new_zern_values.append(
                        x[pos] * zern_n1_m[pos] + y[pos] * zern_n1_n1m[pos] + x[pos] * zern_n1_m1[pos] - zern_n2_m1[pos]
                    )

                    new_gradient_x.append(n * zern_n1_m[pos] + n * zern_n1_m1[pos] + gradient_x_n2_m1[pos])
                    new_gradient_y.append(n * zern_n1_n1m[pos] + gradient_y_n2_m1[pos])

            elif n % 2 == 0 and m == n / 2:
                # Collect Needed Prior Zernike and Gradient Polynomials

                # Index(n - 1, n - 1 - m)
                zern_n1_n1m = zernike_map.get((n - 1, n - 1 - m))[0][1]

                # Index(n - 1, m)
                zern_n1_m = zernike_map.get((n - 1, m))[0][1]

                # Index(n - 1, m - 1)
                zern_n1_m1 = zernike_map.get((n - 1, m - 1))[0][1]

                # Index(n - 2, m - 1)
                zern_n2_m1 = zernike_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_x_n2_m1 = gradient_x_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_y_n2_m1 = gradient_y_map.get((n - 2, m - 1))[0][1]

                # Iterating through every input coordinate point, define the new
                # Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in range(0, num_coords):
                    # Calculate new values
                    new_zern_values.append(
                        2.0 * x[pos] * zern_n1_m[pos] + 2.0 * y[pos] * zern_n1_m1[pos] - zern_n2_m1[pos]
                    )

                    new_gradient_x.append(2.0 * n * zern_n1_m[pos] + gradient_x_n2_m1[pos])
                    new_gradient_y.append(2.0 * n * zern_n1_n1m[pos] + gradient_y_n2_m1[pos])

            else:
                # Collect Needed Prior Zernike and Gradient Polynomials

                # Index(n - 1, n - m)
                zern_n1_nm = zernike_map.get((n - 1, n - m))[0][1]

                # Index(n - 1, n - 1 - m)
                zern_n1_n1m = zernike_map.get((n - 1, n - 1 - m))[0][1]

                # Index(n - 1, m)
                zern_n1_m = zernike_map.get((n - 1, m))[0][1]

                # Index(n - 1, m - 1)
                zern_n1_m1 = zernike_map.get((n - 1, m - 1))[0][1]

                # Index(n - 2, m - 1)
                zern_n2_m1 = zernike_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_x_n2_m1 = gradient_x_map.get((n - 2, m - 1))[0][1]

                # Index(n - 2, m - 1)
                gradient_y_n2_m1 = gradient_y_map.get((n - 2, m - 1))[0][1]

                # Iterating through every input coordinate point, define the new
                # Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in range(0, num_coords):
                    # Calculate new values
                    new_zern_values.append(
                        x[pos] * zern_n1_m[pos]
                        + y[pos] * zern_n1_n1m[pos]
                        + x[pos] * zern_n1_m1[pos]
                        - y[pos] * zern_n1_nm[pos]
                        - zern_n2_m1[pos]
                    )

                    new_gradient_x.append(n * zern_n1_m[pos] + n * zern_n1_m1[pos] + gradient_x_n2_m1[pos])
                    new_gradient_y.append(n * zern_n1_n1m[pos] - n * zern_n1_nm[pos] + gradient_y_n2_m1[pos])

            # Add the new values to our nifty dictionaries

            zernike_map[(n, m)].append((ansi_order, new_zern_values))
            gradient_x_map[(n, m)].append((ansi_order, new_gradient_x))
            gradient_y_map[(n, m)].append((ansi_order, new_gradient_y))
        #  End of inner azimuthal loop
    # End of outer radial order loop

    return zernike_map, gradient_x_map, gradient_y_map